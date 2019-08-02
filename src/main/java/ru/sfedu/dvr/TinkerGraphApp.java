package ru.sfedu.dvr;

import org.apache.tinkerpop.gremlin.structure.Vertex;
import org.apache.tinkerpop.gremlin.tinkergraph.structure.TinkerGraph;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class TinkerGraphApp extends GraphApp {
    private static final Logger LOGGER = LoggerFactory.getLogger(TinkerGraphApp.class);

    /**
     * Constructs a graph app using the given properties.
     * @param fileName location of the properties file
     */
    public TinkerGraphApp(final String fileName) {
        super(fileName);
        this.supportsSchema = true;
        this.supportsTransactions = false;
    }

    @Override
    public void createSchema() {
        LOGGER.info("creating schema");
        final TinkerGraph tinkerGraph = (TinkerGraph) graph;
        // naive check if the schema was previously created
        if (!tinkerGraph.getIndexedKeys(Vertex.class).iterator().hasNext()) {
            tinkerGraph.createIndex("name", Vertex.class);
        }
    }

    public static void main(String[] args) {
        final String fileName = (args != null && args.length > 0) ? args[0] : null;
        final TinkerGraphApp app = new TinkerGraphApp(fileName);
        app.runApp();
    }
}
