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
public class DepmapSangerCrisprConverter extends BioDirectoryConverter
{
    //
    private static final String DATASET_TITLE = "DepMap Sanger CRISPR";
    private static final String DATA_SOURCE_NAME = "DepMap Sanger CRISPR";

    private static final String TAXON_ID = "9606"; // Human Taxon ID

    private static final String CN_CSV_FILE = "gene_effect.csv";

    protected IdResolver rslv;
    private static final Logger LOG = Logger.getLogger(DepmapSangerCrisprConverter.class);

    private Map<String, String> genes = new HashMap<String, String>();
    private Map<String, String> resolvedGenes = new HashMap<String, String>();
    private Map<String, String> unresolvableGenes = new HashMap<String, String>();
    private Map<String, String> cellLines = new HashMap<String, String>();

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
    public DepmapSangerCrisprConverter(ItemWriter writer, Model model) {
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

        processEffect(new FileReader(files.get(CN_CSV_FILE)), geneListArray);

    }

    private Map<String, File> readFilesInDir(File dir) {
        Map<String, File> files = new HashMap<String, File>();
        for (File file : dir.listFiles()) {
            files.put(file.getName(), file);
        }
        return files;
    }

    private void processEffect(Reader reader, ArrayList<String> geneList) throws ObjectStoreException, IOException {
        Iterator<?> lineIter = FormattedTextParser.parseCsvDelimitedReader(reader);
        // header
        String[] firstLine = (String[]) lineIter.next();
        ArrayList<String> genes = new ArrayList<String>();
        for(int i = 1; i < firstLine.length; i++) {
            String formattedGene = firstLine[i].split(" ")[0].trim();
            genes.add(formattedGene);
        }

        //lineIter.next();
        while (lineIter.hasNext()) {
            String[] line = (String[]) lineIter.next();

            String cellLine = line[0];
            for(int i = 1; i < line.length; i++) {
                String effectValue = line[i];
                String theGeneForThisItem = genes.get(i-1);

                if(!geneList.isEmpty()) {
                    if(unresolvableGenes.get(theGeneForThisItem) != null) {
                        continue;
                    }
                    String geneId = getGeneId(theGeneForThisItem);
                    if(geneId == null) {
                        continue;
                    }
                    if(!geneList.contains(geneId)) {
                        continue;
                    }
                }

                Item SangerCrisprItem;

                SangerCrisprItem = createItem("DepMapSangerCrisprGeneEffect");

                SangerCrisprItem.setReference("cellLine", getCellLine(cellLine));

                if(unresolvableGenes.get(theGeneForThisItem) != null) {
                    continue;
                }
                String geneId = getGeneId(theGeneForThisItem);
                if(geneId == null) {
                    continue;
                }

                SangerCrisprItem.setReference("gene", geneId);

                if(!effectValue.isEmpty() && !effectValue.equals("NA")) {
                    LOG.info("SANGER CRISPR 1");
                    SangerCrisprItem.setAttribute("DepmapSangerCrisprGeneEffectValue", effectValue);
                    store(SangerCrisprItem);
                } else {
                    continue;
                }
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
