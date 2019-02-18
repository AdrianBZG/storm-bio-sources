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
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

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

    private Map<String, String> genes = new HashMap<String, String>();
    private Map<String, String> drugs = new HashMap<String, String>();
    private Map<String, String> publications = new HashMap<String, String>();

    private String organismIdentifier; // Not the taxon ID. It references the object that is created into the database.

    /**
     * Constructor
     * @param writer the ItemWriter used to handle the resultant items
     * @param model the Model
     */
    public DgidbDataConverter(ItemWriter writer, Model model) {
        super(writer, model, DATA_SOURCE_NAME, DATASET_TITLE);
    }

    /**
     * 
     *
     * {@inheritDoc}
     */
    public void process(File dataDir) throws Exception {

        Map<String, File> files = readFilesInDir(dataDir);

        organismIdentifier = getOrganism(TAXON_ID);

        processDrugs(new FileReader(files.get(DRUGS_TSV_FILE)));
        processInteractions(new FileReader(files.get(INTERACTIONS_TSV_FILE)));

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

    private void processInteractions(Reader reader) throws ObjectStoreException, IOException {
        Iterator<?> lineIter = FormattedTextParser.parseTabDelimitedReader(reader);
        // Skip header
        lineIter.next();
        while (lineIter.hasNext()) {
            String[] line = (String[]) lineIter.next();

            if(line.length != 10) {
                continue;
            }

            String geneEntrezId = line[2];
            if(geneEntrezId.isEmpty()) {
                continue;
            }
            String chemblId = line[8];

            String interactionType = line[4];
            String pubmedId = line[9];

            if(chemblId.isEmpty()) {
                continue;
            }

            Item interactionItem;
            interactionItem = createItem("DrugInteraction");
            interactionItem.setReference("gene", getGene(geneEntrezId));
            interactionItem.setReference("drug", getDrug(chemblId));
            if(!interactionType.isEmpty()) {
                interactionItem.setAttribute("type", interactionType);
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

    public String getGene(String identifier) {
        String refId = genes.get(identifier);
        if (refId == null) {
            Item gene = createItem("Gene");
            gene.setAttribute("primaryIdentifier", identifier);
            try {
                store(gene);
            } catch (ObjectStoreException e) {
                throw new RuntimeException("failed to store gene with primary identifier: " + identifier, e);
            }
            refId = gene.getIdentifier();
            genes.put(identifier, refId);
        }
        return refId;
    }
}
