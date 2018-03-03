package edu.mit.broad.picard.quality;

import edu.mit.broad.picard.util.DbSnpBitSetUtil;
import net.sf.picard.util.ListMap;
import net.sf.picard.cmdline.*;
import net.sf.picard.filter.DuplicateReadFilter;
import net.sf.picard.filter.NotPrimaryAlignmentFilter;
import net.sf.picard.filter.SamRecordFilter;
import net.sf.picard.io.IoUtil;
import net.sf.picard.metrics.MetricBase;
import net.sf.picard.metrics.MetricsFile;
import net.sf.picard.reference.ReferenceSequenceFileWalker;
import net.sf.picard.util.*;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMReadGroupRecord;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.SequenceUtil;
import net.sf.samtools.util.StringUtil;

import java.io.File;
import java.util.*;

import static java.lang.Math.log10;

/**
 * Class for trying to quantify the CpCG->CpCA error rate.
 */
public class CollectOxoGMetrics extends CommandLineProgram {

	@Usage
	public final String USAGE =
			CommandLineParser.getStandardUsagePreamble(getClass()) +
			"A utility for quantifying the CpCG -> CpCA error rate.";

    @Option(shortName=StandardOptionDefinitions.INPUT_SHORT_NAME,
			doc="Input BAM file for analysis.")
    public File INPUT;

    @Option(shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME,
			doc="Location of output metrics file to write.")
    public File OUTPUT;

    @Option(shortName=StandardOptionDefinitions.REFERENCE_SHORT_NAME,
			doc="Reference sequence to which BAM is aligned.")
    public File REFERENCE_SEQUENCE;

    @Option(doc="An optional list of intervals to restrict analysis to.",
			optional=true)
    public File INTERVALS;
    
    @Option(doc="VCF format dbSNP file, used to exclude regions around known polymorphisms from analysis.",
			optional=true)
    public File DB_SNP;

    @Option(shortName="Q",
			doc="The minimum base quality score for a base to be included in analysis.")
    public int MINIMUM_QUALITY_SCORE = 20;

    @Option(shortName="REFBASE",
            doc="Specific reference base marking SNP artifact on read 2 REF_BASE>ART_BASE.")
    public String REF_BASE = "C";

    @Option(shortName="ARTBASE",
            doc="Specific alternate base marking SNP artifact on read 2 REF_BASE>ART_BASE.")
    public String ART_BASE = "A";

    @Option(shortName="MQ",
			doc="The minimum mapping quality score for a base to be included in analysis.")
    public int MINIMUM_MAPPING_QUALITY = 30;

    @Option(shortName="MIN_INS",
			doc="The minimum insert size for a read to be included in analysis. Set of 0 to allow unpaired reads.")
    public int MINIMUM_INSERT_SIZE = 60;

    @Option(shortName="MAX_INS",
			doc="The maximum insert size for a read to be included in analysis. Set of 0 to allow unpaired reads.")
    public int MAXIMUM_INSERT_SIZE = 600;

    @Option(doc="When available, use original quality scores for filtering.")
    public boolean USE_OQ = true;

    @Option(doc="The number of context bases to include on each side of the assayed G/C base.")
    public int CONTEXT_SIZE = 1;

    @Option(doc="The optional set of sequence contexts to restrict analysis to. If not supplied all contexts are analyzed.")
    public Set<String> CONTEXTS = new HashSet<String>();

    @Option(doc="For debugging purposes: stop after visiting this many sites with at least 1X coverage.")
    public int STOP_AFTER = Integer.MAX_VALUE;

    private final Log log = Log.getInstance(CollectOxoGMetrics.class);
    
    /** Metrics class for outputs. */
    public static final class CpcgMetrics extends MetricBase {
        /** The name of the sample being assayed. */
        public String SAMPLE_ALIAS;
        /** The name of the library being assayed. */
        public String LIBRARY;
        /** The sequence context being reported on. */
        public String CONTEXT;
        /** The total number of sites that had at least one base covering them. */
        public int TOTAL_SITES;
        /** The total number of basecalls observed at all sites. */
        public long TOTAL_BASES;
        /** The number of reference alleles observed as C in read 1 and G in read 2. */
        public long REF_NONOXO_BASES;
        /** The number of reference alleles observed as G in read 1 and C in read 2. */
        public long REF_OXO_BASES;
		/** The total number of reference alleles observed */
        public long REF_TOTAL_BASES;
		/**
		 * The count of observed A basecalls at C reference positions and T basecalls
		 * at G reference bases that are correlated to instrument read number in a way
		 * that rules out oxidation as the cause
		 */
        public long ALT_NONOXO_BASES;
		/**
		 * The count of observed A basecalls at C reference positions and T basecalls
		 * at G reference bases that are correlated to instrument read number in a way
		 * that is consistent with oxidative damage.
		 */
        public long ALT_OXO_BASES;
		/** The oxo error rate, calculated as max(ALT_OXO_BASES - ALT_NONOXO_BASES, 1) / TOTAL_BASES */
        public double OXIDATION_ERROR_RATE;
		/** -10 * log10(OXIDATION_ERROR_RATE) */
        public double OXIDATION_Q;

    }

    /**
     * SAM filter for insert size range.
     */
    static class InsertSizeFilter implements SamRecordFilter {
        final int minInsertSize;
        final int maxInsertSize;

        InsertSizeFilter(final int minInsertSize, final int maxInsertSize) {
            this.minInsertSize = minInsertSize;
            this.maxInsertSize = maxInsertSize;
        }

        @Override public boolean filterOut(final SAMRecord rec) {
            if (rec.getReadPairedFlag()) {
                final int ins = Math.abs(rec.getInferredInsertSize());
                return ins < minInsertSize || ins > maxInsertSize;
            }

            // If the read isn't paired and either min or max is specified filter it out
            return minInsertSize != 0 || maxInsertSize != 0;
        }

        @Override public boolean filterOut(final SAMRecord r1, final SAMRecord r2) {
            return filterOut(r1) || filterOut(r2);
        }
    }

    // Stock main method
    public static void main(final String[] args) {
        new CollectOxoGMetrics().instanceMainWithExit(args);
    }

    @Override
    protected String[] customCommandLineValidation() {
        final int size = 1 + 2*CONTEXT_SIZE;
        final List<String> messages = new ArrayList<String>();
        
        for (final String ctx : CONTEXTS) {
            if (ctx.length() != size) {
                messages.add("Context " + ctx + " is not " + size + " long as implied by CONTEXT_SIZE=" + CONTEXT_SIZE);
            }
            //else if (ctx.charAt(ctx.length() / 2) != REF_BASEB) {
            //    messages.add("Middle base of context sequence " + ctx + " must be C");
            //}
        }


        return messages.isEmpty() ? null : messages.toArray(new String[messages.size()]);
    }

    @Override
    protected int doWork() {
        IoUtil.assertFileIsReadable(INPUT);
        IoUtil.assertFileIsWritable(OUTPUT);
        if (INTERVALS != null) IoUtil.assertFileIsReadable(INTERVALS);
        IoUtil.assertFileIsReadable(REFERENCE_SEQUENCE);
        
        final ReferenceSequenceFileWalker refWalker = new ReferenceSequenceFileWalker(REFERENCE_SEQUENCE);
        final SAMFileReader in = new SAMFileReader(INPUT);
        
        final Set<String> samples   = new HashSet<String>();
        final Set<String> libraries = new HashSet<String>();
        for (final SAMReadGroupRecord rec : in.getFileHeader().getReadGroups()) {
            samples.add(rec.getSample());
            libraries.add(rec.getLibrary());
        }
        
        // Setup the calculators
        final Set<String> contexts = CONTEXTS.isEmpty() ? makeContextStrings(CONTEXT_SIZE) : CONTEXTS;
        final ListMap<String, Calculator> calculators = new ListMap<String, Calculator>();
        for (final String context : contexts) {
            for (final String library : libraries) {
                calculators.add(context, new Calculator(library, context));
            }
        }
        
        // Load up dbSNP if available
        log.info("Loading dbSNP File: " + DB_SNP);
        final DbSnpBitSetUtil dbSnp;
        if (DB_SNP != null) dbSnp = new DbSnpBitSetUtil(DB_SNP, in.getFileHeader().getSequenceDictionary());
        else dbSnp = null;

        // Make an iterator that will filter out funny looking things
        final SamLocusIterator iterator;
        if (INTERVALS != null) {
            final IntervalList intervals = IntervalList.fromFile(INTERVALS);
            intervals.unique();
            iterator = new SamLocusIterator(in, intervals, false);
        }
        else {
            iterator = new SamLocusIterator(in);
        }
        iterator.setEmitUncoveredLoci(false);
        iterator.setMappingQualityScoreCutoff(MINIMUM_MAPPING_QUALITY);
        iterator.setSamFilters(Arrays.asList(
                new NotPrimaryAlignmentFilter(),
                new DuplicateReadFilter(),
                new InsertSizeFilter(MINIMUM_INSERT_SIZE, MAXIMUM_INSERT_SIZE)
        ));

        log.info("Starting iteration.");
        long nextLogTime = 0;
        int sites = 0;

        byte[] REF_BASEa = REF_BASE.getBytes();
        byte REF_BASEb = REF_BASEa[0];
        byte REF_FLIPb =  SequenceUtil.complement(REF_BASEb);

        for (final SamLocusIterator.LocusInfo info : iterator) {
            // Skip dbSNP sites
            final String chrom = info.getSequenceName();
            final int pos      = info.getPosition();
            final int index    = pos-1;
            if (dbSnp != null && dbSnp.isDbSnpSite(chrom, pos)) continue;
            
            // Skip sites at the end of chromosomes 
            final byte[] bases = refWalker.get(info.getSequenceIndex()).getBases();
            if (pos < 3 || pos > bases.length-3) continue;

            // Skip non C-G bases
            final byte base = StringUtil.toUpperCase(bases[index]);
            if (base != REF_BASEb && base != REF_FLIPb) continue;
            
            // Get the context string
            final String context;
            {
                final String tmp = StringUtil.bytesToString(bases, index-CONTEXT_SIZE, 1 + (2*CONTEXT_SIZE)).toUpperCase();
                if (base == REF_BASEb) context = tmp;
                else /* if G */  context = SequenceUtil.reverseComplement(tmp);
            }

            final List<Calculator> calculatorsForContext = calculators.get(context);
            if (calculatorsForContext == null) continue; // happens if we get ambiguous bases in the reference
            for (final Calculator calc : calculatorsForContext) calc.accept(info, base);

            // See if we need to stop
            if (++sites % 100 == 0) {
                final long now = System.currentTimeMillis();
                if (now > nextLogTime) {
                    log.info("Visited " + sites + " sites of interest. Last site: " + chrom + ":" + pos);
                    nextLogTime = now + 60000;
                }
            }
            if (sites >= STOP_AFTER) break;
        }

        final MetricsFile<CpcgMetrics, Integer> file = getMetricsFile();
        for (final List<Calculator> calcs : calculators.values()) {
            for (final Calculator calc : calcs) {
                final CpcgMetrics m = calc.finish();
                m.SAMPLE_ALIAS = StringUtil.join(",", new ArrayList<String>(samples));
                file.addMetric(m);
            }
        }

        file.write(OUTPUT);
        return 0;
    }

    private Set<String> makeContextStrings(final int contextSize) {
        final Set<String> contexts = new HashSet<String>();

        byte[] REF_BASEa = REF_BASE.getBytes();
        byte REF_BASEb = REF_BASEa[0];

        for (final byte[] kmer : generateAllKmers(2*contextSize + 1)) {
            //if (kmer[contextSize] == 'C') {
            if (kmer[contextSize] == REF_BASEb) {
                    contexts.add(StringUtil.bytesToString(kmer));
            }
        }
        
        log.info("Generated " + contexts.size() + " context strings.");
        return contexts;
    }

    /** Generates all possible kmers of length and returns them as byte[]s. */
    private List<byte[]> generateAllKmers(final int length) {
        final List<byte[]> sofar = new LinkedList<byte[]>();
        final byte[] bases = {'A', 'C', 'G', 'T'};

        if (sofar.size() == 0) {
            sofar.add(new byte[length]);
        }

        while (true) {
            final byte[] bs = sofar.remove(0);
            final int indexOfNextBase = findIndexOfNextBase(bs);

            if (indexOfNextBase == -1) {
                sofar.add(bs);
                break;
            }
            else {
                for (final byte b : bases) {
                    final byte[] next = Arrays.copyOf(bs, bs.length);
                    next[indexOfNextBase] = b;
                    sofar.add(next);
                }
            }
        }

        return sofar;
    }

    /** Finds the first non-zero character in the array, or returns -1 if all are non-zero. */
    private int findIndexOfNextBase(final byte[] bs) {
        for (int i=0; i<bs.length; ++i) {
            if (bs[i] == 0) return i;
        }

        return -1;
    }

    /** A little class for counting alleles. */
    private static class Counts {
        int controlA;
        int oxidatedA;
        int controlC;
        int oxidatedC;
        int total() { return controlC + oxidatedC + controlA + oxidatedA; }
    }

    /**
     * Class that calculated CpCG metrics for a specific library.
     */
    private class Calculator {
        private final String library;
        private final String context;

        // Things to be accumulated
        int sites = 0;
        long controlA  = 0;
        long oxidatedA = 0; 
        long controlC  = 0;
        long oxidatedC = 0;

        
        Calculator(final String library, final String context) {
            this.library = library;
            this.context = context;
        }
        
        void accept(final SamLocusIterator.LocusInfo info, final byte refBase) {
            final Counts counts = computeAlleleFraction(info, refBase);

            if (counts.total() > 0) {
                // Things calculated on all sites with coverage
                this.sites++;
                this.controlA  += counts.controlA;
                this.oxidatedA += counts.oxidatedA;
                this.controlC  += counts.controlC;
                this.oxidatedC += counts.oxidatedC;
                ++sites;
            }
        }
        
        CpcgMetrics finish() {
            final CpcgMetrics m = new CpcgMetrics();
            m.LIBRARY              = this.library;
            m.CONTEXT              = this.context;
            m.TOTAL_SITES          = this.sites;
            m.TOTAL_BASES          = this.controlC + + this.oxidatedC + this.controlA + this.oxidatedA;
            m.REF_OXO_BASES        = this.oxidatedC;
            m.REF_NONOXO_BASES     = this.controlC;
            m.REF_TOTAL_BASES      = m.REF_OXO_BASES + m.REF_NONOXO_BASES;
            m.ALT_NONOXO_BASES     = this.controlA;
            m.ALT_OXO_BASES        = this.oxidatedA;

            /**
             * Why do we calculate the oxo error rate using oxidatedA - controlA you ask?  We know that all the
             * bases counted in oxidatedA are consistent with 8-oxo-G damage during shearing, but not all of them
             * will have been caused by this. If we assume that C>A errors caused by other factors will occur randomly
             * with respect to read1/read2, then we should see as many in the 8-oxo-G consistent state as not.  So we
             * assume that controlA is half the story, and remove the other half from oxidatedA.
             */
            m.OXIDATION_ERROR_RATE = Math.max(this.oxidatedA - this.controlA, 1) / (double) m.TOTAL_BASES;
            m.OXIDATION_Q          = -10 * log10(m.OXIDATION_ERROR_RATE);

            return m;
        }

        /**
         *
         */
        private Counts computeAlleleFraction(final SamLocusIterator.LocusInfo info, final byte refBase) {
            final Counts counts = new Counts();

            byte[] REF_BASEa = REF_BASE.getBytes();
            byte[] ART_BASEa = ART_BASE.getBytes();
            byte REF_BASEb = REF_BASEa[0];
            byte ART_BASEb = ART_BASEa[0];
            byte REF_FLIPb =  SequenceUtil.complement(REF_BASEb);
            byte ART_FLIPb =  SequenceUtil.complement(ART_BASEb);

            final byte altBase = (refBase == REF_BASEb) ? (byte) ART_BASEb : (byte) ART_FLIPb;
            
            for (final SamLocusIterator.RecordAndOffset rec : info.getRecordAndPositions()) {
                final byte qual;
                final SAMRecord samrec = rec.getRecord();
                
                if (USE_OQ) {
                    final byte[] oqs = samrec.getOriginalBaseQualities();
                    if (oqs != null) qual = oqs[rec.getOffset()];
                    else qual = rec.getBaseQuality();
                }
                else {
                    qual = rec.getBaseQuality();
                }

                // Skip if below qual, or if library isn't a match
                if (qual < MINIMUM_QUALITY_SCORE) continue;
                if (!this.library.equals(samrec.getReadGroup().getLibrary())) continue;
                
                // Get the read base, and get it in "as read" orientation
                final byte base = rec.getReadBase();
                final byte baseAsRead = samrec.getReadNegativeStrandFlag() ? SequenceUtil.complement(base) : base;
                final int read        = samrec.getReadPairedFlag() && samrec.getSecondOfPairFlag() ? 2 : 1;

                // Figure out how to count the alternative allele. If the damage is caused by oxidation of G
                // during shearing (in non-rnaseq data), then we know that:
                //     G>T observation is always in read 1
                //     C>A observation is always in read 2
                // But if the substitution is from other causes the distribution of A/T across R1/R2 will be
                // random.

                if (base == refBase) {
                    if      (baseAsRead == REF_FLIPb && read == 1) ++counts.oxidatedC;
                    else if (baseAsRead == REF_FLIPb && read == 2) ++counts.controlC;
                    else if (baseAsRead == REF_BASEb && read == 1) ++counts.controlC;
                    else if (baseAsRead == REF_BASEb && read == 2) ++counts.oxidatedC;
                }
                else if (base == altBase)  {
                    if      (baseAsRead == ART_FLIPb && read == 1) ++counts.oxidatedA;
                    else if (baseAsRead == ART_FLIPb && read == 2) ++counts.controlA;
                    else if (baseAsRead == ART_BASEb && read == 1) ++counts.controlA;
                    else if (baseAsRead == ART_BASEb && read == 2) ++counts.oxidatedA;
                }

//                if (base == refBase) {
//                    if      (baseAsRead == 'G' && read == 1) ++counts.oxidatedC;
//                    else if (baseAsRead == 'G' && read == 2) ++counts.controlC;
//                    else if (baseAsRead == 'C' && read == 1) ++counts.controlC;
//                    else if (baseAsRead == 'C' && read == 2) ++counts.oxidatedC;
//                }
//                else if (base == altBase) & (ALT_BASE == 'A') {
//                    if      (baseAsRead == 'T' && read == 1) ++counts.oxidatedA;
//                    else if (baseAsRead == 'T' && read == 2) ++counts.controlA;
//                    else if (baseAsRead == 'A' && read == 1) ++counts.controlA;
//                    else if (baseAsRead == 'A' && read == 2) ++counts.oxidatedA;
//                }
//
            }

            return counts;
        }
    }
}
