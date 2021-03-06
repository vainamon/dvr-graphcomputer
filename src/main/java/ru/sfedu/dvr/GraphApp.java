package ru.sfedu.dvr;

import org.apache.commons.configuration.Configuration;
import org.apache.commons.configuration.ConfigurationException;
import org.apache.commons.configuration.PropertiesConfiguration;
import org.apache.tinkerpop.gremlin.process.computer.ComputerResult;
import org.apache.tinkerpop.gremlin.process.traversal.dsl.graph.GraphTraversalSource;
import org.apache.tinkerpop.gremlin.structure.Graph;
import org.apache.tinkerpop.gremlin.structure.Vertex;
import org.apache.tinkerpop.gremlin.structure.io.IoCore;
import org.apache.tinkerpop.gremlin.structure.io.gryo.GryoMapper;
import org.apache.tinkerpop.gremlin.structure.io.gryo.GryoReader;
import org.apache.tinkerpop.gremlin.structure.io.gryo.GryoWriter;
import org.apache.tinkerpop.gremlin.structure.util.GraphFactory;
import org.apache.tinkerpop.shaded.kryo.Kryo;
import org.javatuples.Triplet;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.awt.*;
import java.io.*;
import java.util.ArrayList;
import java.util.List;

import org.javatuples.Pair;

import javax.vecmath.Color3f;
import javax.vecmath.Point3d;
import javax.vecmath.Point4i;
import javax.vecmath.Vector3d;

import java.awt.image.BufferedImage;
import javax.imageio.ImageIO;

public class GraphApp {
    private static final Logger LOGGER = LoggerFactory.getLogger(GraphApp.class);

    protected String propFileName;
    protected Configuration conf;
    protected Graph graph;
    protected GraphTraversalSource g;
    protected boolean supportsTransactions;
    protected boolean supportsSchema;

    /**
     * Constructs a graph app using the given properties.
     * @param fileName location of the properties file
     */
    public GraphApp(final String fileName) {
        propFileName = fileName;
    }

    /**
     * Opens the graph instance. If the graph instance does not exist, a new
     * graph instance is initialized.
     */
    public GraphTraversalSource openGraph() throws ConfigurationException {
        LOGGER.info("opening graph, properties: " + propFileName);
        conf = new PropertiesConfiguration(propFileName);
        graph = GraphFactory.open(conf);
        g = graph.traversal();
        return g;
    }

    /**
     * Closes the graph instance.
     */
    public void closeGraph() throws Exception {
        LOGGER.info("closing graph");
        try {
            if (g != null) {
                g.close();
            }
            if (graph != null) {
                graph.close();
            }
        } finally {
            g = null;
            graph = null;
        }
    }

    /**
     * Drops the graph instance. The default implementation does nothing.
     */
    public void dropGraph() throws Exception {
    }

    /**
     * Creates the graph schema. The default implementation does nothing.
     */
    public void createSchema() {
    }

    /**
     * Adds the vertices, edges, and properties to the graph.
     */
    public void createElementsAS(byte[] rawData, int xGridSize, int yGridSize, int zGridSize,
                                  double pixelDistance, double sliceDistance,
                                  int volumeDimension) {
        try {
            // naive check if the graph was previously created
            if (g.V(0).has("volume").hasNext()) {
                if (supportsTransactions) {
                    g.tx().rollback();
                }
                return;
            }
            LOGGER.info("creating elements");

            ArrayList<Pair<Vertex, AABB>> vertices = new ArrayList<>();

            for (int i = 0; i < zGridSize - 1; i = i + volumeDimension) {
                for (int j = 0; j < yGridSize - 1; j = j + volumeDimension) {
                    for (int k = 0; k < xGridSize - 1; k = k + volumeDimension) {

                        Point3d min = new Point3d(k * pixelDistance, j * pixelDistance, i * sliceDistance);

                        int maxii = (i + volumeDimension < zGridSize) ? volumeDimension : zGridSize - 1 - i;
                        int maxjj = (j + volumeDimension < yGridSize) ? volumeDimension : yGridSize - 1 - j;
                        int maxkk = (k + volumeDimension < xGridSize) ? volumeDimension : xGridSize - 1 - k;

                        int N = (maxii + 1) * (maxjj + 1) * (maxkk + 1);

                        int [] subVolume = new int[(maxii + 1) * (maxjj + 1) * (maxkk + 1)];

                        Boolean isEmpty = true;
                        Double mean = 0.0;
                        Double variance = 0.0;
                        int maxValue = Integer.MIN_VALUE;

                        for (int ii = 0; ii < maxii + 1; ii++) {
                            for (int jj = 0; jj < maxjj + 1; jj++) {
                                for (int kk = 0; kk < maxkk + 1; kk++) {
                                    int value = rawData[(kk + k) + (jj + j) * xGridSize + (ii + i) * yGridSize * xGridSize] & 0xff;

                                    /* copy value from input volume to node sub volume */
                                    subVolume[kk + jj * (maxkk + 1) + ii * (maxjj + 1) * (maxkk + 1)] = value;

                                    maxValue = Integer.max(value, maxValue);

                                    mean += value / Double.valueOf(N);
                                }
                            }
                        }

                        VolumeProperties volumeProperties = new VolumeProperties();

                        if (maxValue == 0)
                            volumeProperties.isEmpty = true;
                        else {
                            volumeProperties.isEmpty = false;

                            for (int ii = 0; ii < maxii + 1; ii++) {
                                for (int jj = 0; jj < maxjj + 1; jj++) {
                                    for (int kk = 0; kk < maxkk + 1; kk++) {
                                        int value = subVolume[kk + jj * (maxkk + 1) + ii * (maxjj + 1) * (maxkk + 1)];

                                        variance += Math.pow((value - mean) / Double.valueOf(maxValue), 2);
                                    }
                                }
                            }

                            volumeProperties.normalizedDeviation = Math.sqrt(variance / Double.valueOf(N - 1));
                        }

                        Point3d max = new Point3d((k + maxkk) * pixelDistance,
                                (j + maxjj) * pixelDistance, (i + maxii) * sliceDistance);

                        vertices.add(Pair.with(g.addV().property("volume", Triplet.with(new AABB(min, max), subVolume, volumeProperties)).next(),
                                new AABB(min, max)));
                    }
                }
            }

            int zDim = (zGridSize - 1) / volumeDimension + (((zGridSize - 1) % volumeDimension) == 0 ? 0 : 1);
            int yDim = (yGridSize - 1) / volumeDimension + (((yGridSize - 1) % volumeDimension) == 0 ? 0 : 1);
            int xDim = (xGridSize - 1) / volumeDimension + (((xGridSize - 1) % volumeDimension) == 0 ? 0 : 1);

            for (int i = 0; i < zDim; i++) {
                for (int j = 0; j < yDim; j++) {
                    for (int k = 0; k < xDim; k++) {
                        ArrayList<Triplet<Integer,Integer,Integer>> neighbors = new ArrayList<>();

                        neighbors.add(Triplet.with(k - 1, j - 1, i + 1));
                        neighbors.add(Triplet.with(k - 1, j, i + 1));
                        neighbors.add(Triplet.with(k - 1, j + 1, i + 1));
                        neighbors.add(Triplet.with(k - 1, j + 1, i));

                        neighbors.add(Triplet.with(k, j - 1, i + 1));
                        neighbors.add(Triplet.with(k, j, i + 1));
                        neighbors.add(Triplet.with(k, j + 1, i + 1));
                        neighbors.add(Triplet.with(k, j + 1, i));

                        neighbors.add(Triplet.with(k + 1, j, i + 1));
                        neighbors.add(Triplet.with(k + 1, j, i));
                        neighbors.add(Triplet.with(k + 1, j + 1, i + 1));
                        neighbors.add(Triplet.with(k + 1, j + 1, i));
                        neighbors.add(Triplet.with(k + 1, j - 1, i + 1));

                        for (Triplet<Integer,Integer,Integer> neighbor: neighbors) {
                            if ((neighbor.getValue0() >= 0) && (neighbor.getValue0() < xDim)
                                    && ((neighbor.getValue1() >= 0) && (neighbor.getValue1() < yDim))
                                    && (neighbor.getValue2() < zDim)) {

                                final Pair<Vertex, AABB> originV = vertices.get(k + j * xDim + i * xDim * yDim);

                                final Pair<Vertex, AABB> neighborV = vertices.get(neighbor.getValue0()
                                        + neighbor.getValue1() * xDim
                                        + neighbor.getValue2() * xDim * yDim);

                                g.V(originV.getValue0()).as("a")
                                        .V(neighborV.getValue0()).addE("link1").property("aabb", neighborV.getValue1()).from("a")
                                        .next();
                                g.V(neighborV.getValue0()).as("b")
                                        .V(originV.getValue0()).addE("link2").property("aabb", originV.getValue1()).from("b")
                                        .next();
                            }
                        }
                    }
                }
            }

            LOGGER.info("Dim  " + xDim + "; V count " + g.V().count().next() + "; E count " + g.E().count().next());

            if (supportsTransactions) {
                g.tx().commit();
            }

        } catch (Exception e) {
            LOGGER.error(e.getMessage(), e);
            if (supportsTransactions) {
                g.tx().rollback();
            }
        }
    }

    public void createElementsRaw(byte[] rawData, int xGridSize, int yGridSize, int zGridSize,
                               double pixelDistance, double sliceDistance,
                               int volumeDimension) {
        try {
            // naive check if the graph was previously created
            if (g.V(0).has("volume").hasNext()) {
                if (supportsTransactions) {
                    g.tx().rollback();
                }
                return;
            }
            LOGGER.info("creating elements");

            ArrayList<Pair<Vertex, AABB>> vertices = new ArrayList<>();

            for (int i = 0; i < zGridSize - 1; i = i + volumeDimension) {
                for (int j = 0; j < yGridSize - 1; j = j + volumeDimension) {
                    for (int k = 0; k < xGridSize - 1; k = k + volumeDimension) {

                        Point3d min = new Point3d(k * pixelDistance, j * pixelDistance, i * sliceDistance);

                        int maxii = (i + volumeDimension < zGridSize) ? volumeDimension : zGridSize - 1 - i;
                        int maxjj = (j + volumeDimension < yGridSize) ? volumeDimension : yGridSize - 1 - j;
                        int maxkk = (k + volumeDimension < xGridSize) ? volumeDimension : xGridSize - 1 - k;

                        int [] subVolume = new int[(maxii + 1) * (maxjj + 1) * (maxkk + 1)];

                        Boolean isEmpty = true;

                        for (int ii = 0; ii < maxii + 1; ii++) {
                            for (int jj = 0; jj < maxjj + 1; jj++) {
                                for (int kk = 0; kk < maxkk + 1; kk++) {
                                    /* copy value from input volume to node sub volume */
                                    subVolume[kk + jj * (maxkk + 1) + ii * (maxjj + 1) * (maxkk + 1)] =
                                            rawData[(kk + k) + (jj + j) * xGridSize + (ii + i) * yGridSize * xGridSize] & 0xff;

                                    if (subVolume[kk + jj * (maxkk + 1) + ii * (maxjj + 1) * (maxkk + 1)] != 0)
                                        isEmpty = false;
                                }
                            }
                        }

                        Point3d max = new Point3d((k + maxkk) * pixelDistance,
                                (j + maxjj) * pixelDistance, (i + maxii) * sliceDistance);

                        vertices.add(Pair.with(g.addV().property("volume", Triplet.with(new AABB(min, max), subVolume, isEmpty)).next(),
                                new AABB(min, max)));
                    }
                }
            }

            int zDim = (zGridSize - 1) / volumeDimension + (((zGridSize - 1) % volumeDimension) == 0 ? 0 : 1);
            int yDim = (yGridSize - 1) / volumeDimension + (((yGridSize - 1) % volumeDimension) == 0 ? 0 : 1);
            int xDim = (xGridSize - 1) / volumeDimension + (((xGridSize - 1) % volumeDimension) == 0 ? 0 : 1);

            for (int i = 0; i < zDim; i++) {
                for (int j = 0; j < yDim; j++) {
                    for (int k = 0; k < xDim; k++) {
                        ArrayList<Triplet<Integer,Integer,Integer>> neighbors = new ArrayList<>();

                        neighbors.add(Triplet.with(k - 1, j - 1, i + 1));
                        neighbors.add(Triplet.with(k - 1, j, i + 1));
                        neighbors.add(Triplet.with(k - 1, j + 1, i + 1));
                        neighbors.add(Triplet.with(k - 1, j + 1, i));

                        neighbors.add(Triplet.with(k, j - 1, i + 1));
                        neighbors.add(Triplet.with(k, j, i + 1));
                        neighbors.add(Triplet.with(k, j + 1, i + 1));
                        neighbors.add(Triplet.with(k, j + 1, i));

                        neighbors.add(Triplet.with(k + 1, j, i + 1));
                        neighbors.add(Triplet.with(k + 1, j, i));
                        neighbors.add(Triplet.with(k + 1, j + 1, i + 1));
                        neighbors.add(Triplet.with(k + 1, j + 1, i));
                        neighbors.add(Triplet.with(k + 1, j - 1, i + 1));

                        for (Triplet<Integer,Integer,Integer> neighbor: neighbors) {
                            if ((neighbor.getValue0() >= 0) && (neighbor.getValue0() < xDim)
                                    && ((neighbor.getValue1() >= 0) && (neighbor.getValue1() < yDim))
                                    && (neighbor.getValue2() < zDim)) {

                                final Pair<Vertex, AABB> originV = vertices.get(k + j * xDim + i * xDim * yDim);

                                final Pair<Vertex, AABB> neighborV = vertices.get(neighbor.getValue0()
                                        + neighbor.getValue1() * xDim
                                        + neighbor.getValue2() * xDim * yDim);

                                g.V(originV.getValue0()).as("a")
                                        .V(neighborV.getValue0()).addE("link1").property("aabb", neighborV.getValue1()).from("a")
                                        .next();
                                g.V(neighborV.getValue0()).as("b")
                                        .V(originV.getValue0()).addE("link2").property("aabb", originV.getValue1()).from("b")
                                        .next();
                            }
                        }
                    }
                }
            }

            LOGGER.info("Dim  " + xDim + "; V count " + g.V().count().next() + "; E count " + g.E().count().next());

            if (supportsTransactions) {
                g.tx().commit();
            }

        } catch (Exception e) {
            LOGGER.error(e.getMessage(), e);
            if (supportsTransactions) {
                g.tx().rollback();
            }
        }
    }

    public void createElementsEVC(byte[] rawData, int xGridSize, int yGridSize, int zGridSize,
                                  double pixelDistance, double sliceDistance,
                                  int volumeDimension) {
        try {
            // naive check if the graph was previously created
            if (g.V(0).has("volume").hasNext()) {
                if (supportsTransactions) {
                    g.tx().rollback();
                }
                return;
            }
            LOGGER.info("creating elements");

            ArrayList<Pair<Vertex, AABB>> vertices = new ArrayList<>();

            for (int i = 0; i < zGridSize - 1; i = i + volumeDimension) {
                for (int j = 0; j < yGridSize - 1; j = j + volumeDimension) {
                    for (int k = 0; k < xGridSize - 1; k = k + volumeDimension) {

                        ArrayList<EquidistantVolumeCell> volumeCells = new ArrayList<>();

                        Point3d min = new Point3d(k * pixelDistance, j * pixelDistance, i * sliceDistance);

                        int maxii = i, maxjj = j, maxkk = k;

                        for (int ii = i; (ii < i + volumeDimension) && (ii < zGridSize - 1); ii++) {
                            for (int jj = j; (jj < j + volumeDimension) && (jj < yGridSize - 1) ; jj++) {
                                for (int kk = k; (kk < k + volumeDimension) && (kk < xGridSize - 1) ; kk++) {

                                    Point3d origin = new Point3d(kk * pixelDistance, jj * pixelDistance, ii * sliceDistance);

                                    // V000 -> V100 -> V010 -> V110
                                    Point4i frontIntensities = new Point4i(
                                            rawData[kk + jj * xGridSize + ii * yGridSize * xGridSize] & 0xff,
                                            rawData[(kk + 1) + jj * xGridSize + ii * yGridSize * xGridSize] & 0xff,
                                            rawData[kk + (jj + 1) * xGridSize + ii * yGridSize * xGridSize] & 0xff,
                                            rawData[(kk + 1) + (jj + 1) * xGridSize + ii * yGridSize * xGridSize] & 0xff
                                    );

                                    // V001 -> V101 -> V011 -> V111
                                    Point4i backIntensities = new Point4i(
                                            rawData[kk + jj * xGridSize + (ii + 1) * yGridSize * xGridSize] & 0xff,
                                            rawData[(kk + 1) + jj * xGridSize + (ii + 1) * yGridSize * xGridSize] & 0xff,
                                            rawData[kk + (jj + 1) * xGridSize + (ii + 1) * yGridSize * xGridSize] & 0xff,
                                            rawData[(kk + 1) + (jj + 1) * xGridSize + (ii + 1) * yGridSize * xGridSize] & 0xff
                                    );

                                    EquidistantVolumeCell volumeCell = new EquidistantVolumeCell(origin,
                                            frontIntensities, backIntensities);

                                    volumeCells.add(volumeCell);

                                    maxkk = kk > maxkk ? kk : maxkk;
                                }

                                maxjj = jj > maxjj ? jj : maxjj;
                            }

                            maxii = ii > maxii ? ii : maxii;
                        }

                        Point3d max = new Point3d((maxkk + 1) * pixelDistance,
                                (maxjj + 1) * pixelDistance, (maxii + 1) * sliceDistance);

                        vertices.add(Pair.with(g.addV().property("volume", Pair.with(new AABB(min, max), volumeCells)).next(),
                                new AABB(min, max)));
                    }
                }
            }

            int zDim = (zGridSize - 1) / volumeDimension + (((zGridSize - 1) % volumeDimension) == 0 ? 0 : 1);
            int yDim = (yGridSize - 1) / volumeDimension + (((yGridSize - 1) % volumeDimension) == 0 ? 0 : 1);
            int xDim = (xGridSize - 1) / volumeDimension + (((xGridSize - 1) % volumeDimension) == 0 ? 0 : 1);

            for (int i = 0; i < zDim; i++) {
                for (int j = 0; j < yDim; j++) {
                    for (int k = 0; k < xDim; k++) {
                        ArrayList<Triplet<Integer,Integer,Integer>> neighbors = new ArrayList<>();

                        neighbors.add(Triplet.with(k - 1, j - 1, i + 1));
                        neighbors.add(Triplet.with(k - 1, j, i + 1));
                        neighbors.add(Triplet.with(k - 1, j + 1, i + 1));
                        neighbors.add(Triplet.with(k - 1, j + 1, i));

                        neighbors.add(Triplet.with(k, j - 1, i + 1));
                        neighbors.add(Triplet.with(k, j, i + 1));
                        neighbors.add(Triplet.with(k, j + 1, i + 1));
                        neighbors.add(Triplet.with(k, j + 1, i));

                        neighbors.add(Triplet.with(k + 1, j, i + 1));
                        neighbors.add(Triplet.with(k + 1, j, i));
                        neighbors.add(Triplet.with(k + 1, j + 1, i + 1));
                        neighbors.add(Triplet.with(k + 1, j + 1, i));
                        neighbors.add(Triplet.with(k + 1, j - 1, i + 1));

                        for (Triplet<Integer,Integer,Integer> neighbor: neighbors) {
                            if ((neighbor.getValue0() >= 0) && (neighbor.getValue0() < xDim)
                                    && ((neighbor.getValue1() >= 0) && (neighbor.getValue1() < yDim))
                                    && (neighbor.getValue2() < zDim)) {

                                final Pair<Vertex, AABB> originV = vertices.get(k + j * xDim + i * xDim * yDim);

                                final Pair<Vertex, AABB> neighborV = vertices.get(neighbor.getValue0()
                                        + neighbor.getValue1() * xDim
                                        + neighbor.getValue2() * xDim * yDim);

                                g.V(originV.getValue0()).as("a")
                                        .V(neighborV.getValue0()).addE("link1").property("aabb", neighborV.getValue1()).from("a")
                                        .next();
                                g.V(neighborV.getValue0()).as("b")
                                        .V(originV.getValue0()).addE("link2").property("aabb", originV.getValue1()).from("b")
                                        .next();
                            }
                        }
                    }
                }
            }

            LOGGER.info("Dim  " + xDim + "; V count " + g.V().count().next() + "; E count " + g.E().count().next());

            if (supportsTransactions) {
                g.tx().commit();
            }

        } catch (Exception e) {
            LOGGER.error(e.getMessage(), e);
            if (supportsTransactions) {
                g.tx().rollback();
            }
        }
    }

    public void render(int width, int height, double pixelDistance, double sliceDistance, double samplingStep, AABB volumeBox) {
        Camera camera = new Camera();

        Point3d cameraOrigin = new Point3d(0,0,-100);
        Point3d target = new Point3d(0,0,0);
        Vector3d up = new Vector3d(0,1,0);

        camera.set(cameraOrigin, target, up);
        camera.slide(-32,-32, -50);
        //camera.slide(-105,-128, 20);
        //camera.slide(-345,-340, 320); /* zeiss */

        Point3d P1 = new Point3d(4,3,-5);
        Point3d P2 = new Point3d(-4,3,-5);
        Point3d P3 = new Point3d(4,-3,-5);
        Point3d P4 = new Point3d(-4,-3,-5);

        Point3d origin = new Point3d(0, 0, 0);

        origin.set(camera.transform(origin));

        P1.set(camera.transform(P1));
        P2.set(camera.transform(P2));
        P3.set(camera.transform(P3));
        P4.set(camera.transform(P4));

        Vector3d DX = new Vector3d(P2);

        DX.sub(P1);
        DX.scale(1.0 / (width - 1));

        Vector3d DY = new Vector3d(P3);

        DY.sub(P1);
        DY.scale(1.0 / (height - 1));

        LOGGER.info("Scene AABB: " + volumeBox);

        try {
            if (g == null) {
                return;
            }

            DVRVertexProgramAS dvrvp = DVRVertexProgramAS.build()
                    .height(height)
                    .width(width)
                    .P1(P1)
                    .DX(DX)
                    .DY(DY)
                    .origin(origin)
                    //.samplingStep(samplingStep)
                    .minSamplingStep(samplingStep / 5000)
                    .maxSamplingStep(samplingStep)
                    .p(2)
                    .pixelDistance(pixelDistance)
                    .sliceDistance(sliceDistance)
                    .sceneBox(volumeBox)
                    .create();

//            DVRVertexProgram dvrvp = DVRVertexProgram.build()
//                    .height(height)
//                    .width(width)
//                    .P1(P1)
//                    .DX(DX)
//                    .DY(DY)
//                    .origin(origin)
//                    .samplingStep(samplingStep)
//                    .pixelDistance(pixelDistance)
//                    .sliceDistance(sliceDistance)
//                    .sceneBox(volumeBox)
//                    .create();
            ComputerResult result = graph.compute().program(dvrvp).submit().get();

            LOGGER.info("DVR runtime = " + result.memory().getRuntime() + "ms; iteration: " + result.memory().getIteration());

            if (result.memory().exists(DVRVertexProgramAS.IMAGE_COLORS_OUT_PARAM_STRING)) {
                List<Pair<Pair<Integer, Integer>, Color3f>> colors = result.memory().get(DVRVertexProgramAS.IMAGE_COLORS_OUT_PARAM_STRING);

                final BufferedImage image =
                        new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);

                for (int i = 0; i < width; i++) {
                    for (int j = 0; j < height; j++) {
                        Color c = new Color(100, 125, 120);
                        image.setRGB(i, j, c.getRGB());
                    }
                }

                for (Pair<Pair<Integer, Integer>, Color3f> resultPair : colors) {
                    int r = Math.min(Math.round(resultPair.getValue1().getX() * 255), 255);
                    int g = Math.min(Math.round(resultPair.getValue1().getY() * 255), 255);
                    int b = Math.min(Math.round(resultPair.getValue1().getZ() * 255), 255);

                    Color resultColor = new Color(r, g, b);
                    image.setRGB(resultPair.getValue0().getValue0(), resultPair.getValue0().getValue1(), resultColor.getRGB());
                }

                File outputfile = new File("saved.png");
                ImageIO.write(image, "png", outputfile);
            }

            if (supportsTransactions) {
                g.tx().commit();
            }
        } catch (Exception e) {
            LOGGER.error(e.getMessage(), e);
            if (supportsTransactions) {
                g.tx().rollback();
            }
        }
    }

    public void runApp() {
        try {
            // open and initialize the graph
            openGraph();

            // define the schema before loading data
            if (supportsSchema) {
                createSchema();
            }

            int xGridSize = 64, yGridSize = 64, zGridSize = 64;
//            int xGridSize = 512, yGridSize = 256, zGridSize = 256;
//            int xGridSize = 256, yGridSize = 256, zGridSize = 256;
//            int xGridSize = 680, yGridSize = 680, zGridSize = 680;
            double pixelDistance = 1.0, sliceDistance = 1.0;

            ClassLoader classLoader = Thread.currentThread().getContextClassLoader();

            FileInputStream fis = new FileInputStream(new File(classLoader.getResource("neghip_64x64x64_uint8.raw").getFile()));
//            FileInputStream fis = new FileInputStream(new File(classLoader.getResource("golfball0_0-512x256x256.raw").getFile()));
//            FileInputStream fis = new FileInputStream(new File(classLoader.getResource("bonsai_256x256x256_uint8.raw").getFile()));
//            FileInputStream fis = new FileInputStream(new File(classLoader.getResource("zeiss_680x680x680_uint8.raw").getFile()));
            byte [] buffer = new byte[fis.available()];
            fis.read(buffer);
            fis.close();

            // build the graph structure
//            createElementsEVC(buffer, xGridSize, yGridSize, zGridSize, pixelDistance, sliceDistance, 4);
//              createElementsRaw(buffer, xGridSize, yGridSize, zGridSize, pixelDistance, sliceDistance, 32);
            createElementsAS(buffer, xGridSize, yGridSize, zGridSize, pixelDistance, sliceDistance, 2);

//            File file = new File("zeiss-graph-128.kryo");
//            OutputStream fos = new FileOutputStream(file);
//            GryoMapper mapper = GryoMapper.build()
//                    .addCustom(ru.sfedu.dvr.AABB.class)
//                    .addCustom(javax.vecmath.Point3d.class)
//                    .addCustom(ru.sfedu.dvr.VolumeProperties.class).create();
//            GryoWriter writer = GryoWriter.build().mapper(mapper).create();
//            writer.writeGraph(fos, graph);

//            File file = new File("golfball-graph.kryo");
//            InputStream fis = new FileInputStream(file);
//            GryoMapper mapper = GryoMapper.build()
//                    .addCustom(ru.sfedu.dvr.AABB.class)
//                    .addCustom(javax.vecmath.Point3d.class)
//                    .addCustom(ru.sfedu.dvr.VolumeProperties.class).create();
//            GryoReader reader = GryoReader.build().mapper(mapper).create();
//            reader.readGraph(fis, graph);

            Point3d min = new Point3d(0, 0, 0);

            Point3d max = new Point3d((xGridSize - 1) * pixelDistance,
                    (yGridSize - 1) * pixelDistance,
                    (zGridSize - 1) * sliceDistance);

            AABB volumeBox = new AABB(min, max);

            render(1024, 768, pixelDistance, sliceDistance, 0.5, volumeBox);

            // close the graph
            closeGraph();
        } catch (Exception e) {
            LOGGER.error(e.getMessage(), e);
        }
    }
}
