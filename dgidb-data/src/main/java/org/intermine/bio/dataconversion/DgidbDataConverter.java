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
public class DgidbDataConverter extends BioDirectoryConverter
{
    //
    private static final String DATASET_TITLE = "DGIdb Dataset";
    private static final String DATA_SOURCE_NAME = "DGIdb";

    private static final String TAXON_ID = "9606"; // Human Taxon ID

    private static final String DRUGS_TSV_FILE = "drugs.tsv";
    private static final String INTERACTIONS_TSV_FILE = "interactions.tsv";

    protected IdResolver rslv;

    private Map<String, String> genes = new HashMap<String, String>();
    private Map<String, String> drugs = new HashMap<String, String>();

    private static final Logger LOG = Logger.getLogger(DgidbDataConverter.class);

    private Map<String, String> publications = new HashMap<String, String>();

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
    public DgidbDataConverter(ItemWriter writer, Model model) {
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

        processDrugs(new FileReader(files.get(DRUGS_TSV_FILE)));
        processInteractions(new FileReader(files.get(INTERACTIONS_TSV_FILE)), geneListArray);

    }

    private Map<String, File> readFilesInDir(File dir) {
        Map<String, File> files = new HashMap<String, File>();
        for (File file : dir.listFiles()) {
            files.put(file.getName(), file);
        }
        return files;
    }

    private void processDrugs(Reader reader) throws ObjectStoreException, IOException {
        Iterator<?> lineIter = FormattedTextParser.parseTabDelimitedReader(reader);
        // Skip header
        lineIter.next();
        while (lineIter.hasNext()) {
            String[] line = (String[]) lineIter.next();

            // Check for missing values in the line
            if(line.length != 4) {
                continue;
            }

            String drugName = line[1];
            String chemblId = line[2];
            String source = line[3];

            // Check for missing values in the line
            if(chemblId.isEmpty()) {
                continue;
            }

            Item drugItem;
            drugItem = createItem("Drug");
            drugItem.setAttribute("primaryIdentifier", chemblId);
            drugItem.setAttribute("name", drugName);
            drugItem.setAttribute("source", source);

            store(drugItem);
            drugs.put(chemblId, drugItem.getIdentifier());

            // Look Chembl API for more information?
        }
    }

    private void processInteractions(Reader reader, ArrayList<String> geneList) throws ObjectStoreException, IOException {
        Iterator<?> lineIter = FormattedTextParser.parseTabDelimitedReader(reader);
        // Skip header
        lineIter.next();
        while (lineIter.hasNext()) {
            String[] line = (String[]) lineIter.next();

            if(line.length != 10) {
                continue;
            }

            String geneSymbol = line[0];
            if(geneSymbol.isEmpty()) {
                continue;
            }

            if(!geneList.isEmpty()) {
                String resolvedGene = getGeneIdentifier(geneSymbol);
                if(!geneList.contains(resolvedGene)) {
                    continue;
                }
            }

            String chemblId = line[8];

            String interactionType = line[4];
            String pubmedId = line[9];

            if(chemblId.isEmpty()) {
                continue;
            }

            String geneId = getGeneId(geneSymbol);

            if (StringUtils.isEmpty(geneId)) {
                continue;
            }

            Item interactionItem;
            interactionItem = createItem("DrugInteraction");
            interactionItem.setReference("gene", geneId);
            interactionItem.setReference("drug", getDrug(chemblId));
            if(!StringUtils.isEmpty(interactionType)) {
                interactionItem.setAttribute("type", interactionType);
            } else {
                continue;
            }

            if(!pubmedId.isEmpty()) {
                interactionItem.setReference("publication", getPublication(pubmedId));
            }

            store(interactionItem);
        }
    }

    public String getPublication(String identifier) {
        String refId = publications.get(identifier);
        if (refId == null) {
            Item publication = createItem("Publication");
            publication.setAttribute("pubMedId", identifier);
            try {
                store(publication);
            } catch (ObjectStoreException e) {
                throw new RuntimeException("failed to store publication with pubMedId: " + identifier, e);
            }
            refId = publication.getIdentifier();
            publications.put(identifier, refId);
        }
        return refId;
    }

    public String getDrug(String identifier) {
        String refId = drugs.get(identifier);
        if (refId == null) {
            Item drug = createItem("Drug");
            drug.setAttribute("primaryIdentifier", identifier);
            try {
                store(drug);
            } catch (ObjectStoreException e) {
                throw new RuntimeException("failed to store drug with primary identifier: " + identifier, e);
            }
            refId = drug.getIdentifier();
            drugs.put(identifier, refId);
        }
        return refId;
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
}
