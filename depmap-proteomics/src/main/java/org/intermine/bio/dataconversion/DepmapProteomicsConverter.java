package org.intermine.bio.dataconversion;

/*
 * Copyright (C) 2002-2019 FlyMine
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
public class DepmapProteomicsConverter extends BioDirectoryConverter
{
    private static final String DATASET_TITLE = "DepMap CCLE Proteomics Data";
    private static final String DATA_SOURCE_NAME = "DepMap CCLE";

    private static final String TAXON_ID = "9606"; // Human Taxon ID

    private static final String CSV_FILE = "protein_quant_current_normalized.csv";

    private Map<String, String> genes = new HashMap<String, String>();
    private Map<String, String> cellLines = new HashMap<String, String>();
    private Map<String, String> proteins = new HashMap<String, String>();

    protected IdResolver rslv;
    private static final Logger LOG = Logger.getLogger(DepmapProteomicsConverter.class);

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

    public DepmapProteomicsConverter(ItemWriter writer, Model model) {
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

        processData(new FileReader(files.get(CSV_FILE)), geneListArray);
    }

    private Map<String, File> readFilesInDir(File dir) {
        Map<String, File> files = new HashMap<String, File>();
        for (File file : dir.listFiles()) {
            files.put(file.getName(), file);
        }
        return files;
    }

    private void processData(Reader reader, ArrayList<String> geneList) throws ObjectStoreException, IOException {
        Iterator<?> lineIter = FormattedTextParser.parseCsvDelimitedReader(reader);

        // header
        String[] firstLine = (String[]) lineIter.next();
        ArrayList<String> cellLinesList = new ArrayList<String>();
        for(int i = 48; i < firstLine.length; i++) {
            String formattedCL = firstLine[i].split("_Ten")[0].trim();
            cellLinesList.add(formattedCL);
        }

        while (lineIter.hasNext()) {
            String[] line = (String[]) lineIter.next();

            String geneSymbol = line[1];

            if(!geneList.isEmpty()) {
                String resolvedGene = getGeneIdentifier(geneSymbol);
                if(!geneList.contains(resolvedGene)) {
                    continue;
                }
            }

            String proteinID = line[0];
            String uniprotID = line[4];

            // For each of the cell lines in the row
            for(int i = 48; i < line.length; i++) {
                String theCLForThisItem = cellLinesList.get(i-48);
                String proteinQuant = line[i];

                Item integratedItem = createItem("DepMapProteomics");

                if(!theCLForThisItem.isEmpty()) {
                    integratedItem.setReference("cellLine", getCellLine(theCLForThisItem));
                } else {
                    continue;
                }

                if(!geneSymbol.isEmpty()) {
                    String geneId = getGeneId(geneSymbol);

                    if (StringUtils.isEmpty(geneId)) {
                        continue;
                    }

                    integratedItem.setReference("gene", geneId);
                } else {
                    continue;
                }

                if(!proteinID.isEmpty()) {
                    integratedItem.setAttribute("ProteinID", proteinID);
                } else {
                    continue;
                }

                if(!uniprotID.isEmpty()) {
                    integratedItem.setAttribute("ProteinUniprotID", uniprotID);
                    //integratedItem.setReference("ProteinUniprotID", getProtein(uniprotID));
                } else {
                    continue;
                }

                if(!proteinQuant.isEmpty() && !proteinQuant.equals("NA")) {
                    integratedItem.setAttribute("ProteinQuantitation", proteinQuant);
                } else {
                    continue;
                }

                store(integratedItem);
                
            }
        }
    }

    private String getGeneId(String primaryIdentifier) throws ObjectStoreException {
        String resolvedIdentifier = resolveGene(primaryIdentifier);
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

    private String resolveGene(String identifier) {
        String id = identifier;

        if (rslv != null && rslv.hasTaxon(TAXON_ID)) {
            int resCount = rslv.countResolutions(TAXON_ID, identifier);
            if (resCount != 1) {
                LOG.info("RESOLVER: failed to resolve gene to one identifier, ignoring gene: "
                        + identifier + " count: " + resCount + " Human identifier: "
                        + rslv.resolveId(TAXON_ID, identifier));
                return null;
            }
            id = rslv.resolveId(TAXON_ID, identifier).iterator().next();
        }
        return id;
    }

    public String getProtein(String identifier) {
        String refId = proteins.get(identifier);
        if (refId == null) {
            Item pr = createItem("Protein");
            pr.setAttribute("uniprotName", identifier);
            try {
                store(pr);
            } catch (ObjectStoreException e) {
                throw new RuntimeException("failed to store Protein with uniprotName: " + identifier, e);
            }
            refId = pr.getIdentifier();
            proteins.put(identifier, refId);
        }
        return refId;
    }

    public String getCellLine(String identifier) {
        String refId = cellLines.get(identifier);
        if (refId == null) {
            Item cl = createItem("CellLine");
            cl.setAttribute("CCLEname", identifier);
            try {
                store(cl);
            } catch (ObjectStoreException e) {
                throw new RuntimeException("failed to store cell line with CCLEname: " + identifier, e);
            }
            refId = cl.getIdentifier();
            cellLines.put(identifier, refId);
        }
        return refId;
    }
}
