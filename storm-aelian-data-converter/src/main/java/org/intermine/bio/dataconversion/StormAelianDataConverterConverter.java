package org.intermine.bio.dataconversion;

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
public class StormAelianDataConverterConverter extends BioDirectoryConverter
{
    private static final String DATASET_TITLE = "StormAelianData";
    private static final String DATA_SOURCE_NAME = "StormAelianData";

    private static final String TAXON_ID = "9606"; // Human Taxon ID

    private static final String CSV_FILE = "WTA_markers_after_guide_enrichment_no_threshold.csv";

    private Map<String, String> genes = new HashMap<String, String>();
    private Map<String, String> resolvedGenes = new HashMap<String, String>();
    private Map<String, String> unresolvableGenes = new HashMap<String, String>();

    protected IdResolver rslv;
    private static final Logger LOG = Logger.getLogger(StormAelianDataConverterConverter.class);

    private String organismIdentifier; // Not the taxon ID. It references the object that is created into the database.
    //

    public StormAelianDataConverterConverter(ItemWriter writer, Model model) {
        super(writer, model, DATA_SOURCE_NAME, DATASET_TITLE);
        if (rslv == null) {
            rslv = IdResolverService.getIdResolverByOrganism(TAXON_ID);
        }
    }

    public void process(File dataDir) throws Exception {
        Map<String, File> files = readFilesInDir(dataDir);

        organismIdentifier = getOrganism(TAXON_ID);

        processData(new FileReader(files.get(CSV_FILE)));
    }

    private Map<String, File> readFilesInDir(File dir) {
        Map<String, File> files = new HashMap<String, File>();
        for (File file : dir.listFiles()) {
            files.put(file.getName(), file);
        }
        return files;
    }

    private void processData(Reader reader) throws ObjectStoreException, IOException {
        Iterator<?> lineIter = FormattedTextParser.parseCsvDelimitedReader(reader);

        // header
        String[] firstLine = (String[]) lineIter.next();

        while (lineIter.hasNext()) {
            String[] line = (String[]) lineIter.next();

            for(int i = 0; i < line.length; i++) {
                String marker = line[0];
                String p_val = line[1];
                String avg_log2FC = line[2];
                String p_val_adj = line[5];            
                String ident = line[6];

                Item integratedItem = createItem("StormAelianData");

                String markerId = getGeneId(marker);
                if(markerId == null) {
                    continue;
                }

                integratedItem.setReference("marker", markerId);

                if(!p_val.isEmpty()) {
                    integratedItem.setAttribute("p_val", p_val);
                } else {
                    continue;
                }

                if(!avg_log2FC.isEmpty()) {
                    integratedItem.setAttribute("avg_log2FC", avg_log2FC);
                } else {
                    continue;
                }

                if(!p_val_adj.isEmpty()) {
                    integratedItem.setAttribute("p_val_adj", p_val_adj);
                } else {
                    continue;
                }

                String identId = getGeneId(ident);
                if(identId == null) {
                    continue;
                }

                integratedItem.setReference("ident", identId);

                store(integratedItem);
                
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

    private boolean isDouble(String str) {
        try {
            // check if it can be parsed as any double
            double x = Double.parseDouble(str);
            // check if the double can be converted without loss to an int
            if (x == (int) x)
                // if yes, this is an int, thus return false
                return false;
            // otherwise, this cannot be converted to an int (e.g. "1.2")
            return true;
            // short version: return x != (int) x;
        }
        catch(NumberFormatException e) {
            return false;
        }
    }
}