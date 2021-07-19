package org.intermine.bio.dataconversion;

/*
 * Copyright (C) 2002-2018 FlyMine
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  See the LICENSE file for more
 * information or http://www.gnu.org/copyleft/lesser.html.
 *
 */

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.Reader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

import org.apache.commons.lang.StringUtils;
import org.apache.log4j.Logger;
import org.intermine.dataconversion.ItemWriter;
import org.intermine.metadata.Model;
import org.intermine.objectstore.ObjectStoreException;
import org.intermine.util.FormattedTextParser;
import org.intermine.xml.full.Item;


/**
 * 
 * @author
 */
public class DepmapCcleMutationsConverter extends BioDirectoryConverter
{
    //
    private static final String DATASET_TITLE = "DepMap CCLE Mutations Data";
    private static final String DATA_SOURCE_NAME = "DepMap Public 20Q3";

    private static final String TAXON_ID = "9606"; // Human Taxon ID

    private static final String MUTATIONS_CSV_FILE = "CCLE_mutations.csv";

    private Map<String, String> genes = new HashMap<String, String>();
    private Map<String, String> resolvedGenes = new HashMap<String, String>();
    private Map<String, String> unresolvableGenes = new HashMap<String, String>();
    private Map<String, String> cellLines = new HashMap<String, String>();

    protected IdResolver rslv;
    private static final Logger LOG = Logger.getLogger(DepmapCcleMutationsConverter.class);

    private String organismIdentifier; // Not the taxon ID. It references the object that is created into the database.

    // Methods to integrate the data only for a list of genes
    private static final String GENE_LIST_FILE = "/data/storm/targets/storm_targets_symbols.csv";
    private ArrayList<String> processGeneList(String geneListFile) throws Exception {
        File geneListF = new File(geneListFile);

        Iterator<?> lineIter = FormattedTextParser.parseCsvDelimitedReader(new FileReader(geneListF));
        ArrayList<String> geneListArray = new ArrayList<String>();

        while (lineIter.hasNext()) {
            String[] line = (String[]) lineIter.next();
            String gene = line[0];
            if(StringUtils.isEmpty(gene)) {
                continue;
            }

            if(unresolvableGenes.get(gene) != null) {
                continue;
            }
            String geneId = getGeneId(gene);
            if(geneId == null) {
                continue;
            }
            geneListArray.add(geneId);
        }

        return geneListArray;
    }
    //

    /**
     * Constructor
     * @param writer the ItemWriter used to handle the resultant items
     * @param model the Model
     */
    public DepmapCcleMutationsConverter(ItemWriter writer, Model model) {
        super(writer, model, DATA_SOURCE_NAME, DATASET_TITLE);
        if (rslv == null) {
            rslv = IdResolverService.getIdResolverByOrganism(TAXON_ID);
        }
    }

    public void process(File dataDir) throws Exception {
        Map<String, File> files = readFilesInDir(dataDir);

        organismIdentifier = getOrganism(TAXON_ID);

        ArrayList<String> geneListArray = new ArrayList<String>();
        if(!StringUtils.isEmpty(GENE_LIST_FILE)) {
            geneListArray = processGeneList(GENE_LIST_FILE);
        }

        processMutationsData(new FileReader(files.get(MUTATIONS_CSV_FILE)), geneListArray);
    }

    private Map<String, File> readFilesInDir(File dir) {
        Map<String, File> files = new HashMap<String, File>();
        for (File file : dir.listFiles()) {
            files.put(file.getName(), file);
        }
        return files;
    }

    private void processMutationsData(Reader reader, ArrayList<String> geneList) throws ObjectStoreException, IOException {
        Iterator<?> lineIter = FormattedTextParser.parseTabDelimitedReader(reader);

        // Skip header
        lineIter.next();

        while (lineIter.hasNext()) {
            String[] line = (String[]) lineIter.next();

            Item mutationsItem;

            mutationsItem = createItem("DepMapMutations");

            /*
                -	Hugo Symbol 0
                -	Chromosome 3
                -	Start position 4
                -	End position 5
                -	Strand 6
                -	Variant classification 7
                -	Variant type 8
                -	Genome change 13
                -	Annotation transcript 14
                -	isDeleterious 19
                -	isTCGAhotspot 20
                -	TCGAhsCnt 21
                -	isCOSMIChotspot 22
                -	COSMIChsCnt 23
                -	Variant annotation 32
                -	DepMap ID 33
            */

            String hugoSymbol = line[0];

            if(!geneList.isEmpty()) {
                if(unresolvableGenes.get(hugoSymbol) != null) {
                    continue;
                }
                String geneId = getGeneId(hugoSymbol);
                if(geneId == null) {
                    continue;
                }
                if(!geneList.contains(geneId)) {
                    continue;
                }
            }

            String chromosome = line[3];
            String startPosition = line[4];
            String endPosition = line[5];
            String strand = line[6];
            String variantClassification = line[7];
            String variantType = line[8];
            String genomeChange = line[13];
            String annotationTranscript = line[14];
            String isDeleterious = line[19];
            String isTCGAhotspot = line[20];
            String TCGAhsCnt = line[21];
            String isCOSMIChotspot = line[22];
            String COSMIChsCnt = line[23];
            String variantAnnotation = line[32];
            String DepMapID = line[33];



            if(!DepMapID.isEmpty()) {
                mutationsItem.setReference("cellLine", getCellLine(DepMapID));
            } else {
                continue;
            }

            if(!hugoSymbol.isEmpty()) {
                if(unresolvableGenes.get(hugoSymbol) != null) {
                    continue;
                }
                String geneId = getGeneId(hugoSymbol);
                if(geneId == null) {
                    continue;
                }

                mutationsItem.setReference("gene", geneId);
            } else {
                continue;
            }

            if(!chromosome.isEmpty()) {
                mutationsItem.setAttribute("Chromosome", chromosome);
            } else {
                mutationsItem.setAttribute("Chromosome", "Not specified");
            }

            if(!startPosition.isEmpty()) {
                mutationsItem.setAttribute("Start", startPosition);
            } else {
                mutationsItem.setAttribute("Start", "Not specified");
            }

            if(!endPosition.isEmpty()) {
                mutationsItem.setAttribute("End", endPosition);
            } else {
                mutationsItem.setAttribute("End", "Not specified");
            }

            if(!strand.isEmpty()) {
                mutationsItem.setAttribute("Strand", strand);
            } else {
                mutationsItem.setAttribute("Strand", "Not specified");
            }

            if(!variantClassification.isEmpty()) {
                mutationsItem.setAttribute("VariantClassification", variantClassification);
            } else {
                mutationsItem.setAttribute("VariantClassification", "Not specified");
            }

            if(!variantType.isEmpty()) {
                mutationsItem.setAttribute("VariantType", variantType);
            } else {
                mutationsItem.setAttribute("VariantType", "Not specified");
            }

            if(!genomeChange.isEmpty()) {
                mutationsItem.setAttribute("GenomeChange", genomeChange);
            } else {
                mutationsItem.setAttribute("GenomeChange", "Not specified");
            }

            if(!annotationTranscript.isEmpty()) {
                mutationsItem.setAttribute("AnnotationTranscript", annotationTranscript);
            } else {
                mutationsItem.setAttribute("AnnotationTranscript", "Not specified");
            }

            if(!isDeleterious.isEmpty()) {
                mutationsItem.setAttribute("isDeleterious", isDeleterious);
            } else {
                mutationsItem.setAttribute("isDeleterious", "Not specified");
            }

            if(!isTCGAhotspot.isEmpty()) {
                mutationsItem.setAttribute("isTCGAhotspot", isTCGAhotspot);
            } else {
                mutationsItem.setAttribute("isTCGAhotspot", "Not specified");
            }

            if(!TCGAhsCnt.isEmpty()) {
                mutationsItem.setAttribute("TCGAhsCnt", TCGAhsCnt);
            } else {
                mutationsItem.setAttribute("TCGAhsCnt", "Not specified");
            }

            if(!isCOSMIChotspot.isEmpty()) {
                mutationsItem.setAttribute("isCOSMIChotspot", isCOSMIChotspot);
            } else {
                mutationsItem.setAttribute("isCOSMIChotspot", "Not specified");
            }

            if(!COSMIChsCnt.isEmpty()) {
                mutationsItem.setAttribute("COSMIChsCnt", COSMIChsCnt);
            } else {
                mutationsItem.setAttribute("COSMIChsCnt", "Not specified");
            }

            if(!variantAnnotation.isEmpty()) {
                mutationsItem.setAttribute("VariantAnnotation", variantAnnotation);
            } else {
                mutationsItem.setAttribute("VariantAnnotation", "Not specified");
            }

            store(mutationsItem);
        }
    }

    private String getGeneId(String identifier) throws ObjectStoreException {
        String geneId = null;
        try {
            String resolvedIdentifier = resolveGene(identifier);
            if(resolvedIdentifier != null) {
                geneId = genes.get(resolvedIdentifier);
                if (geneId == null) {
                    Item gene = createItem("Gene");
                    gene.setAttribute("primaryIdentifier", resolvedIdentifier);
                    store(gene);
                    geneId = gene.getIdentifier();
                    genes.put(resolvedIdentifier, geneId);
                }
                return geneId;
            } else {
                return resolvedIdentifier;
            }
        } catch (Exception e) {
            LOG.info("getGeneId: failed to resolve gene: " + identifier);
            return null;
        }
    }

    private String resolveGene(String identifier) {
        String id = null;

        if(resolvedGenes.get(identifier) != null) {
            id = resolvedGenes.get(identifier);
        } else {
            if (rslv != null && rslv.hasTaxon(TAXON_ID)) {
                int resCount = rslv.countResolutions(TAXON_ID, identifier);
                if (resCount != 1) {
                    unresolvableGenes.put(identifier, identifier);
                    return null;
                }
                id = rslv.resolveId(TAXON_ID, identifier).iterator().next();
                resolvedGenes.put(identifier, id);
            }
        }
        return id;
    }

    public String getCellLine(String identifier) {
        String refId = cellLines.get(identifier);
        if (refId == null) {
            Item cl = createItem("CellLine");
            cl.setAttribute("DepMapID", identifier);
            try {
                store(cl);
            } catch (ObjectStoreException e) {
                throw new RuntimeException("failed to store cell line with DepMapID: " + identifier, e);
            }
            refId = cl.getIdentifier();
            cellLines.put(identifier, refId);
        }
        return refId;
    }
}
