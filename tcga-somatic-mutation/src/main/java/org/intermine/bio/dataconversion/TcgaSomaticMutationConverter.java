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
public class TcgaSomaticMutationConverter extends BioDirectoryConverter
{
    //
    private static final String DATASET_TITLE = "TCGA Mutation Data";
    private static final String DATA_SOURCE_NAME = "TCGA Gene-level non-silent mutation PANCAN";

    private static final String TAXON_ID = "9606"; // Human Taxon ID

    private static final String MUTATION_TSV_FILE = "mc3.v0.2.8.PUBLIC.nonsilentGene.xena";

    private Map<String, String> genes = new HashMap<String, String>();
    private Map<String, String> resolvedGenes = new HashMap<String, String>();
    private Map<String, String> unresolvableGenes = new HashMap<String, String>();
    private Map<String, String> samples = new HashMap<String, String>();

    protected IdResolver rslv;
    private static final Logger LOG = Logger.getLogger(TcgaSomaticMutationConverter.class);

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

            String resolvedGeneIdentifier = getGeneIdentifier(gene);
            if(resolvedGeneIdentifier != null) {
                geneListArray.add(resolvedGeneIdentifier);
            }
        }

        return geneListArray;
    }

    private String getGeneIdentifier(String geneSymbol) throws ObjectStoreException {
        String resolvedIdentifier = resolveGene(geneSymbol);
        if (StringUtils.isEmpty(resolvedIdentifier)) {
            return null;
        }
        String geneId = genes.get(resolvedIdentifier);
        if (geneId == null) {
            Item gene = createItem("Gene");
            gene.setAttribute("primaryIdentifier", resolvedIdentifier);
            //gene.setAttribute("symbol", primaryIdentifier);
            //gene.setReference("organism", getOrganism(TAXON_ID));
            store(gene);
            geneId = gene.getIdentifier();
            genes.put(resolvedIdentifier, geneId);
        }
        return geneId;
    }
    //

    /**
     * Constructor
     * @param writer the ItemWriter used to handle the resultant items
     * @param model the Model
     */
    public TcgaSomaticMutationConverter(ItemWriter writer, Model model) {
        super(writer, model, DATA_SOURCE_NAME, DATASET_TITLE);
        if (rslv == null) {
            rslv = IdResolverService.getIdResolverByOrganism(TAXON_ID);
        }
    }

    /**
     * 
     *
     * {@inheritDoc}
     */
    public void process(File dataDir) throws Exception {
        Map<String, File> files = readFilesInDir(dataDir);

        organismIdentifier = getOrganism(TAXON_ID);

        ArrayList<String> geneListArray = new ArrayList<String>();
        if(!StringUtils.isEmpty(GENE_LIST_FILE)) {
            geneListArray = processGeneList(GENE_LIST_FILE);
        }

        processMutationData(new FileReader(files.get(MUTATION_TSV_FILE)), geneListArray);

    }

    private Map<String, File> readFilesInDir(File dir) {
        Map<String, File> files = new HashMap<String, File>();
        for (File file : dir.listFiles()) {
            files.put(file.getName(), file);
        }
        return files;
    }

    private void processMutationData(Reader reader, ArrayList<String> geneList) throws ObjectStoreException, IOException {
        Iterator<?> lineIter = FormattedTextParser.parseTabDelimitedReader(reader);
        // header
        String[] firstLine = (String[]) lineIter.next();
        ArrayList<String> samples = new ArrayList<String>();
        for(int i = 1; i < firstLine.length; i++) {
            String formattedSample = firstLine[i].split(" ")[0].trim();
            samples.add(formattedSample);
        }

        //lineIter.next();
        while (lineIter.hasNext()) {
            String[] line = (String[]) lineIter.next();

            String gene = line[0];

            if(!geneList.isEmpty()) {
                String resolvedGene = getGeneIdentifier(gene);
                if(!geneList.contains(resolvedGene)) {
                    continue;
                }
            }

            for(int i = 1; i < line.length; i++) {
                String mutationValue = line[i];
                String theSampleForThisItem = samples.get(i-1);
                Item MutationItem;

                MutationItem = createItem("TCGAMutation");

                if(!gene.isEmpty()) {
                    if(unresolvableGenes.get(gene) != null) {
                        continue;
                    }
                    String geneId = getGeneId(gene);
                    if(geneId == null) {
                        continue;
                    }

                    MutationItem.setReference("gene", geneId);
                } else {
                    continue;
                }

                if(!theSampleForThisItem.isEmpty()) {
                    MutationItem.setReference("sample", getSample(theSampleForThisItem));
                } else {
                    continue;
                }

                if(!mutationValue.isEmpty() && StringUtils.isNumeric(mutationValue)) {
                    MutationItem.setAttribute("TcgaSomaticMutationValue", mutationValue);
                } else {
                    continue;
                }

                store(MutationItem);
                //cellLines.put(cellLine, CopyNumberItem.getIdentifier());
            }
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

    public String getSample(String identifier) {
        String refId = samples.get(identifier);
        if (refId == null) {
            Item cl = createItem("TCGASample");
            cl.setAttribute("SampleID", identifier);
            try {
                store(cl);
            } catch (ObjectStoreException e) {
                throw new RuntimeException("failed to store sample with ID: " + identifier, e);
            }
            refId = cl.getIdentifier();
            samples.put(identifier, refId);
        }
        return refId;
    }
}
