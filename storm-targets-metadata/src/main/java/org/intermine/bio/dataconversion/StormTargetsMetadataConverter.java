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
public class StormTargetsMetadataConverter extends BioDirectoryConverter
{
    //
    private static final String DATASET_TITLE = "STORM Targets Metadata";
    private static final String DATA_SOURCE_NAME = "Metadata describing the list of STORM targets";

    private static final String TAXON_ID = "9606"; // Human Taxon ID

    private static final String TARGETS_FILE = "storm_targets.csv";
    private static final String CATEGORIES_FILE = "storm_targets_categories.csv";

    private Map<String, String> genes = new HashMap<String, String>();
    private Map<String, String> typeMap = new HashMap<String, String>();

    protected IdResolver rslv;
    private static final Logger LOG = Logger.getLogger(StormTargetsMetadataConverter.class);

    private String organismIdentifier;

    public StormTargetsMetadataConverter(ItemWriter writer, Model model) {
        super(writer, model, DATA_SOURCE_NAME, DATASET_TITLE);
        if (rslv == null) {
            rslv = IdResolverService.getIdResolverByOrganism(TAXON_ID);
        }
    }

    public void process(File dataDir) throws Exception {
        Map<String, File> files = readFilesInDir(dataDir);

        organismIdentifier = getOrganism(TAXON_ID);

        processTargetsCategories(new FileReader(files.get(CATEGORIES_FILE)));
        processTargetsMetadata(new FileReader(files.get(TARGETS_FILE)));
    }

    private Map<String, File> readFilesInDir(File dir) {
        Map<String, File> files = new HashMap<String, File>();
        for (File file : dir.listFiles()) {
            files.put(file.getName(), file);
        }
        return files;
    }

    private void processTargetsCategories(Reader reader) throws ObjectStoreException, IOException {
        Iterator<?> lineIter = FormattedTextParser.parseCsvDelimitedReader(reader);
        // header
        String[] firstLine = (String[]) lineIter.next();

        //lineIter.next();
        while (lineIter.hasNext()) {
            String[] line = (String[]) lineIter.next();

            String typeAbbreviation = line[0];
            String typeName = line[1];

            if (StringUtils.isEmpty(typeAbbreviation) || StringUtils.isEmpty(typeName)) {
                continue;
            }

            String typeNameFromMap = typeMap.get(typeAbbreviation);
            if (typeNameFromMap == null) {
                typeMap.put(typeAbbreviation, typeName);
            }

        }
    }

    private void processTargetsMetadata(Reader reader) throws ObjectStoreException, IOException {
        Iterator<?> lineIter = FormattedTextParser.parseCsvDelimitedReader(reader);
        // header
        String[] firstLine = (String[]) lineIter.next();

        //lineIter.next();
        while (lineIter.hasNext()) {
            String[] line = (String[]) lineIter.next();

            String gene = line[0];
            String modification = line[8];
            String typeAbbreviation = line[9];

            if (StringUtils.isEmpty(typeAbbreviation) || typeMap.get(typeAbbreviation) == null) {
                continue;
            }

            String typeName = typeMap.get(typeAbbreviation);

            String notes = line[10];

            Item IntegratedItem;

            IntegratedItem = createItem("STORMTarget");

            if(!gene.isEmpty()) {
                String geneId = getGeneId(gene);

                if (StringUtils.isEmpty(geneId)) {
                    continue;
                }

                IntegratedItem.setReference("gene", geneId);
            } else {
                continue;
            }

            if(!StringUtils.isEmpty(modification)) {
                IntegratedItem.setAttribute("Modification", modification);
            }
            IntegratedItem.setAttribute("TypeAbbreviation", typeAbbreviation);
            IntegratedItem.setAttribute("TypeName", typeName);
            if(!StringUtils.isEmpty(notes)) {
                IntegratedItem.setAttribute("Notes", notes);
            }

            store(IntegratedItem);
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
}
