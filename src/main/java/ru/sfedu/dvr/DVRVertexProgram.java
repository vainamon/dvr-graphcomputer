/*
 * Licensed to the Apache Software Foundation (ASF) under one
 * or more contributor license agreements.  See the NOTICE file
 * distributed with this work for additional information
 * regarding copyright ownership.  The ASF licenses this file
 * to you under the Apache License, Version 2.0 (the
 * "License"); you may not use this file except in compliance
 * with the License.  You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing,
 * software distributed under the License is distributed on an
 * "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
 * KIND, either express or implied.  See the License for the
 * specific language governing permissions and limitations
 * under the License.
 */

package ru.sfedu.dvr;

import org.apache.commons.configuration.Configuration;
import org.apache.tinkerpop.gremlin.process.computer.GraphComputer;
import org.apache.tinkerpop.gremlin.process.computer.Memory;
import org.apache.tinkerpop.gremlin.process.computer.MemoryComputeKey;
import org.apache.tinkerpop.gremlin.process.computer.MessageScope;
import org.apache.tinkerpop.gremlin.process.computer.Messenger;
import org.apache.tinkerpop.gremlin.process.computer.VertexComputeKey;
import org.apache.tinkerpop.gremlin.process.computer.VertexProgram;
import org.apache.tinkerpop.gremlin.process.computer.util.AbstractVertexProgramBuilder;
import org.apache.tinkerpop.gremlin.process.traversal.Traversal;
import org.apache.tinkerpop.gremlin.process.traversal.Operator;
import org.apache.tinkerpop.gremlin.process.traversal.dsl.graph.__;
import org.apache.tinkerpop.gremlin.process.traversal.util.PureTraversal;
import org.apache.tinkerpop.gremlin.structure.Edge;
import org.apache.tinkerpop.gremlin.structure.Graph;
import org.apache.tinkerpop.gremlin.structure.Vertex;
import org.apache.tinkerpop.gremlin.structure.util.StringFactory;
import org.javatuples.Pair;
import org.javatuples.Triplet;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import javax.vecmath.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.function.Function;
import java.util.Iterator;

public class DVRVertexProgram implements VertexProgram<DVRVertexProgram.PhotonMessage> {

    private static final Logger LOGGER = LoggerFactory.getLogger(DVRVertexProgram.class);

    @SuppressWarnings("WeakerAccess")
    public static final String IMAGE_COLORS_OUT_PARAM_STRING = "ru.sfedu.test.DVRVertexProgram.imageColors";

    private static final String P1_IN_PARAM_STRING = "ru.sfedu.test.DVRVertexProgram.P1";
    private static final String DX_IN_PARAM_STRING = "ru.sfedu.test.DVRVertexProgram.DX";
    private static final String DY_IN_PARAM_STRING = "ru.sfedu.test.DVRVertexProgram.DY";
    private static final String ORIGIN_IN_PARAM_STRING = "ru.sfedu.test.DVRVertexProgram.origin";
    private static final String HEIGHT_IN_PARAM_STRING = "ru.sfedu.test.DVRVertexProgram.height";
    private static final String WIDTH_IN_PARAM_STRING = "ru.sfedu.test.DVRVertexProgram.width";
    private static final String SAMPLING_STEP_IN_PARAM_STRING = "ru.sfedu.test.DVRVertexProgram.samplingStep";
    private static final String PIXEL_DISTANCE_IN_PARAM_STRING = "ru.sfedu.test.DVRVertexProgram.pixelDistance";
    private static final String SLICE_DISTANCE_IN_PARAM_STRING = "ru.sfedu.test.DVRVertexProgram.sliceDistance";
    private static final String SCENE_BOX_IN_PARAM_STRING = "ru.sfedu.test.DVRVertexProgram.sceneBox";

    private static final String VOTE_TO_HALT = "ru.sfedu.test.DVRVertexProgram.voteToHalt";

    public static final PureTraversal<Vertex, Edge> DEFAULT_EDGE_TRAVERSAL = new PureTraversal<>(__.outE().asAdmin());

    private PureTraversal<Vertex, Edge> edgeTraversal = DEFAULT_EDGE_TRAVERSAL.clone();

    private Point3d P1;
    private Vector3d DX;
    private Vector3d DY;
    private Point3d origin;
    private Integer height;
    private Integer width;
    private Number pixelDistance;
    private Number sliceDistance;
    private Number samplingStep;
    private AABB sceneBox;

    private final long time = System.currentTimeMillis();

    private final double EPSILON = 1e-6;

    private final Color3f ambientColor = new Color3f(0.0f, 0.0f, 0.0f);

    private float kA = 4.0f;

    private static final Set<VertexComputeKey> VERTEX_COMPUTE_KEYS = new HashSet<>(Arrays.asList(
            // VertexComputeKey.of()
    ));

    private final Set<MemoryComputeKey> memoryComputeKeys = new HashSet<>(Arrays.asList(
            MemoryComputeKey.of(VOTE_TO_HALT, Operator.and, false, true)
    ));

    private DVRVertexProgram() {

    }

    @Override
    public void loadState(final Graph graph, final Configuration configuration) {
        if (configuration.containsKey(P1_IN_PARAM_STRING))
            this.P1 = (Point3d) configuration.getProperty(P1_IN_PARAM_STRING);
        else
            this.P1 = new Point3d();

        if (configuration.containsKey(DX_IN_PARAM_STRING))
            this.DX = (Vector3d) configuration.getProperty(DX_IN_PARAM_STRING);
        else
            this.DX = new Vector3d();

        if (configuration.containsKey(DY_IN_PARAM_STRING))
            this.DY = (Vector3d) configuration.getProperty(DY_IN_PARAM_STRING);
        else
            this.DY = new Vector3d();

        if (configuration.containsKey(ORIGIN_IN_PARAM_STRING))
            this.origin = (Point3d) configuration.getProperty(ORIGIN_IN_PARAM_STRING);
        else
            this.origin = new Point3d();

        if (configuration.containsKey(SCENE_BOX_IN_PARAM_STRING))
            this.sceneBox = (AABB) configuration.getProperty(SCENE_BOX_IN_PARAM_STRING);
        else
            this.sceneBox = new AABB(new Point3d(), new Point3d());

        if (configuration.containsKey(SAMPLING_STEP_IN_PARAM_STRING))
            this.samplingStep = (Number) configuration.getProperty(SAMPLING_STEP_IN_PARAM_STRING);
        else
            this.samplingStep = 0.1;

        if (configuration.containsKey(PIXEL_DISTANCE_IN_PARAM_STRING))
            this.pixelDistance = (Number) configuration.getProperty(PIXEL_DISTANCE_IN_PARAM_STRING);
        else
            this.pixelDistance = 1.0;

        if (configuration.containsKey(SLICE_DISTANCE_IN_PARAM_STRING))
            this.sliceDistance = (Number) configuration.getProperty(SLICE_DISTANCE_IN_PARAM_STRING);
        else
            this.sliceDistance = 1.0;


        this.height = configuration.getInteger(HEIGHT_IN_PARAM_STRING, 480);
        this.width = configuration.getInteger(WIDTH_IN_PARAM_STRING, 640);

        this.memoryComputeKeys.add(MemoryComputeKey.of(IMAGE_COLORS_OUT_PARAM_STRING, Operator.addAll, true, false));
    }

    @Override
    public void storeState(final Configuration configuration) {
        VertexProgram.super.storeState(configuration);

        configuration.setProperty(HEIGHT_IN_PARAM_STRING, this.height);
        configuration.setProperty(WIDTH_IN_PARAM_STRING, this.width);

        if (this.P1 != null)
            configuration.setProperty(P1_IN_PARAM_STRING, P1);
        if (this.DX != null)
            configuration.setProperty(DX_IN_PARAM_STRING, DX);
        if (this.DY != null)
            configuration.setProperty(DY_IN_PARAM_STRING, DY);
        if (this.origin != null)
            configuration.setProperty(ORIGIN_IN_PARAM_STRING, origin);
        if (this.sceneBox != null)
            configuration.setProperty(SCENE_BOX_IN_PARAM_STRING, sceneBox);
        if (this.samplingStep != null)
            configuration.setProperty(SAMPLING_STEP_IN_PARAM_STRING, samplingStep);
        if (this.pixelDistance != null)
            configuration.setProperty(PIXEL_DISTANCE_IN_PARAM_STRING, pixelDistance);
        if (this.sliceDistance != null)
            configuration.setProperty(SLICE_DISTANCE_IN_PARAM_STRING, sliceDistance);
    }

    @Override
    public Set<VertexComputeKey> getVertexComputeKeys() {
        return VERTEX_COMPUTE_KEYS;
    }

    @Override
    public Set<MemoryComputeKey> getMemoryComputeKeys() {
        return memoryComputeKeys;
    }

    @Override
    public Set<MessageScope> getMessageScopes(final Memory memory) {
        return Collections.emptySet();
    }

    @Override
    public VertexProgram<PhotonMessage> clone() {
        try {
            final DVRVertexProgram clone = (DVRVertexProgram) super.clone();

            if (null != this.edgeTraversal)
                clone.edgeTraversal = this.edgeTraversal.clone();

            return clone;
        } catch (final CloneNotSupportedException e) {
            throw new IllegalStateException(e.getMessage(), e);
        }
    }

    @Override
    public GraphComputer.ResultGraph getPreferredResultGraph() {
        return GraphComputer.ResultGraph.ORIGINAL;
    }

    @Override
    public GraphComputer.Persist getPreferredPersist() {
        return GraphComputer.Persist.NOTHING;
    }

    @Override
    public void setup(final Memory memory) {
        memory.set(VOTE_TO_HALT, true);
    }

    @Override
    public void execute(final Vertex vertex, final Messenger<PhotonMessage> messenger, final Memory memory) {

        boolean voteToHalt = true;

        Triplet<AABB, int[], Boolean> volume =
                (Triplet<AABB, int[], Boolean>) vertex.property("volume").value();

        if (memory.isInitialIteration()) {
            // является ли ячейка объема, хранимая в вершине, граничной,
            // т.е. содержащей точки, принадлежащие плоскостям ограничивающего всю сцену объема
            if (checkIfBoundaryBox(volume.getValue0())) {
                for (int i = 0; i < height * width; i++) {
                    int y = i / width;
                    int x = i % width;

                    Point3d viewPlanePoint = new Point3d();

                    viewPlanePoint.scaleAdd(y, DY, P1);
                    viewPlanePoint.scaleAdd(x, DX, viewPlanePoint);

                    Vector3d rayDir = new Vector3d(viewPlanePoint);

                    rayDir.sub(origin);

                    rayDir.normalize();

                    Ray r = new Ray(origin, rayDir);

                    Point3d q = new Point3d();

                    // определяем пересечение луча с ячейкой объема, хранимой в вершине
                    boolean intersect = volume.getValue0().intersectRay(r, q);

                    if (intersect) {
                        Point3d p = new Point3d();

                        // проверка принадлежности точки пересечения граничным плоскостям объема сцены
                        sceneBox.intersectRay(r, p);

                        if (q.distance(p) <= EPSILON) {
                            PhotonMessage photon = new PhotonMessage("Msg" + i, r, i);
                            photon.setAccColor(ambientColor);

                            traceRay(photon, volume);

                            if (!photon.isInsignificant()) {
                                if (forwardPhoton(photon, vertex, messenger))
                                    voteToHalt = false;
                                else
                                    collectPhoton(photon, memory);
                            } else
                                collectPhoton(photon, memory);
                        }
                    }
                }
            }
        } else {
            final Iterator<PhotonMessage> photonsIterator = messenger.receiveMessages();

            while (photonsIterator.hasNext()) {
                final PhotonMessage nextPhoton = photonsIterator.next();

                traceRay(nextPhoton, volume);

                if (!nextPhoton.isInsignificant()) {
                    if (forwardPhoton(nextPhoton, vertex, messenger))
                        voteToHalt = false;
                    else
                        collectPhoton(nextPhoton, memory);
                } else
                    collectPhoton(nextPhoton, memory);
            }
        }

        memory.add(VOTE_TO_HALT, voteToHalt);
    }

    @Override
    public boolean terminate(final Memory memory) {
        final boolean voteToHalt = memory.get(VOTE_TO_HALT);

        memory.set(VOTE_TO_HALT, true);

        if (voteToHalt)
            return true;
        else
            return false;
    }

    @Override
    public String toString() {

        final List<String> options = new ArrayList<>();
        final Function<String, String> shortName = name -> name.substring(name.lastIndexOf(".") + 1);

        if (!this.edgeTraversal.equals(DEFAULT_EDGE_TRAVERSAL)) {
            options.add(shortName.apply("EDGE_TRAVERSAL") + "=" + this.edgeTraversal.get());
        }

        return StringFactory.vertexProgramString(this, String.join(", ", options));
    }

    //////////////////////////////

    private boolean checkIfBoundaryBox(final AABB volumeBox) {
        if ((Math.abs(sceneBox.min().getX() - volumeBox.min().getX()) <= EPSILON) ||
                (Math.abs(sceneBox.min().getY() - volumeBox.min().getY()) <= EPSILON) ||
                (Math.abs(sceneBox.min().getZ() - volumeBox.min().getZ()) <= EPSILON) ||
                (Math.abs(sceneBox.max().getX() - volumeBox.max().getX()) <= EPSILON) ||
                (Math.abs(sceneBox.max().getY() - volumeBox.max().getY()) <= EPSILON) ||
                (Math.abs(sceneBox.max().getZ() - volumeBox.max().getZ()) <= EPSILON))
            return true;

        return false;
    }

    private void traceRay(PhotonMessage photon, final Triplet<AABB, int[], Boolean> subVolume) {
        Point3d q = new Point3d();

        Double currentSamplePointT;

        Ray ray = photon.ray();

        Double lastSamplePointT = photon.lastSamplePointT();

        subVolume.getValue0().intersectRay(ray, q);

        if (lastSamplePointT == null)
            lastSamplePointT = new Double(ray.getTNear());

        currentSamplePointT = lastSamplePointT;

        Double exitPointT = new Double(ray.getTFar());

        Boolean isEmpty = subVolume.getValue2();

        if (isEmpty) {
            currentSamplePointT += (Math.floor((exitPointT - currentSamplePointT) / samplingStep.doubleValue()) + 1) * samplingStep.doubleValue();
        } else {
            Point3d currentSamplePoint = ray.getPointByT(currentSamplePointT);

            while (exitPointT > currentSamplePointT) {
                int cellMinIndices[] = {
                        (int) Math.floor(Math.abs(currentSamplePoint.getX() - subVolume.getValue0().min().getX()) / pixelDistance.doubleValue()),
                        (int) Math.floor(Math.abs(currentSamplePoint.getY() - subVolume.getValue0().min().getY()) / pixelDistance.doubleValue()),
                        (int) Math.floor(Math.abs(currentSamplePoint.getZ() - subVolume.getValue0().min().getZ()) / sliceDistance.doubleValue())};

                int subVolumeSizes[] = {
                        (int) Math.round(Math.abs(subVolume.getValue0().max().getX() - subVolume.getValue0().min().getX()) / pixelDistance.doubleValue()) + 1,
                        (int) Math.round(Math.abs(subVolume.getValue0().max().getY() - subVolume.getValue0().min().getY()) / pixelDistance.doubleValue()) + 1,
                        (int) Math.round(Math.abs(subVolume.getValue0().max().getZ() - subVolume.getValue0().min().getZ()) / sliceDistance.doubleValue()) + 1};

                for (int i = 0; i < 3; i++)
                    if (cellMinIndices[i] >= subVolumeSizes[i] - 1)
                        cellMinIndices[i] = subVolumeSizes[i] - 2;

                Point3d cellMin = new Point3d(
                        subVolume.getValue0().min().getX() + cellMinIndices[0] * pixelDistance.doubleValue(),
                        subVolume.getValue0().min().getY() + cellMinIndices[1] * pixelDistance.doubleValue(),
                        subVolume.getValue0().min().getZ() + cellMinIndices[2] * sliceDistance.doubleValue()
                );

                Point3d cellMax = new Point3d(
                        subVolume.getValue0().min().getX() + (cellMinIndices[0] + 1) * pixelDistance.doubleValue(),
                        subVolume.getValue0().min().getY() + (cellMinIndices[1] + 1) * pixelDistance.doubleValue(),
                        subVolume.getValue0().min().getZ() + (cellMinIndices[2] + 1) * sliceDistance.doubleValue()
                );

                int subVolumeIntensities[] = subVolume.getValue1();

                // V000 -> V100 -> V010 -> V110
                Point4i fIs = new Point4i(
                        subVolumeIntensities[cellMinIndices[0] + cellMinIndices[1] * subVolumeSizes[0] + cellMinIndices[2] * subVolumeSizes[1] * subVolumeSizes[0]],
                        subVolumeIntensities[cellMinIndices[0] + 1 + cellMinIndices[1] * subVolumeSizes[0] + cellMinIndices[2] * subVolumeSizes[1] * subVolumeSizes[0]],
                        subVolumeIntensities[cellMinIndices[0] + (cellMinIndices[1] + 1) * subVolumeSizes[0] + cellMinIndices[2] * subVolumeSizes[1] * subVolumeSizes[0]],
                        subVolumeIntensities[cellMinIndices[0] + 1 + (cellMinIndices[1] + 1) * subVolumeSizes[0] + cellMinIndices[2] * subVolumeSizes[1] * subVolumeSizes[0]]
                );
                // V001 -> V101 -> V011 -> V111
                Point4i bIs = new Point4i(
                        subVolumeIntensities[cellMinIndices[0] + cellMinIndices[1] * subVolumeSizes[0] + (cellMinIndices[2] + 1) * subVolumeSizes[1] * subVolumeSizes[0]],
                        subVolumeIntensities[cellMinIndices[0] + 1 + cellMinIndices[1] * subVolumeSizes[0] + (cellMinIndices[2] + 1) * subVolumeSizes[1] * subVolumeSizes[0]],
                        subVolumeIntensities[cellMinIndices[0] + (cellMinIndices[1] + 1) * subVolumeSizes[0] + (cellMinIndices[2] + 1) * subVolumeSizes[1] * subVolumeSizes[0]],
                        subVolumeIntensities[cellMinIndices[0] + 1 + (cellMinIndices[1] + 1) * subVolumeSizes[0] + (cellMinIndices[2] + 1) * subVolumeSizes[1] * subVolumeSizes[0]]
                );

                AABB cellAABB = new AABB(cellMin, cellMax);

                cellAABB.intersectRay(ray, q);

                while (ray.getTFar() >= currentSamplePointT) {
                    if (cellAABB.sqDistanceToPoint(currentSamplePoint) < EPSILON) {
                        double xWeight = Math.abs((currentSamplePoint.getX() - cellMin.getX()) / pixelDistance.doubleValue());
                        double yWeight = Math.abs((currentSamplePoint.getY() - cellMin.getY()) / pixelDistance.doubleValue());
                        double zWeight = Math.abs((currentSamplePoint.getZ() - cellMin.getZ()) / sliceDistance.doubleValue());

                        Double interpolatedValue = ((fIs.getX() * (1 - xWeight) + fIs.getY() * xWeight) * (1 - yWeight)
                                + (fIs.getZ() * (1 - xWeight) + fIs.getW() * xWeight) * yWeight) * (1 - zWeight)
                                + ((bIs.getX() * (1 - xWeight) + bIs.getY() * xWeight) * (1 - yWeight)
                                + (bIs.getZ() * (1 - xWeight) + bIs.getW() * xWeight) * yWeight) * zWeight;

                        photon.accumulateCT(phongShading(ctf(interpolatedValue)), otf(interpolatedValue));

                        if (photon.isInsignificant())
                            break;
                    }

                    currentSamplePointT += samplingStep.doubleValue();

                    currentSamplePoint = ray.getPointByT(currentSamplePointT);
                }

                if (photon.isInsignificant())
                    break;
            }
        }

        photon.setLastSamplePointT(currentSamplePointT);
        photon.setExitPointT(exitPointT);
    }

    private Color3f ctf(Double volumeValue) {
        if (volumeValue < EPSILON)
            return new Color3f(0.0f, 0.0f, 0.0f);

        return new Color3f(0.005f, 0.001f, 0.01f);
    }

    private double otf(Double volumeValue) {
        return Math.min(volumeValue / 100.0, 1.0);
    }

    private Color3f phongShading(Color3f classifiedColor) {
        Color3f shadedColor = new Color3f(classifiedColor);

        shadedColor.scale(kA);
        shadedColor.clampMax(1.0f);

        return shadedColor;
    }

    private boolean forwardPhoton(PhotonMessage photon, final Vertex vertex, final Messenger<PhotonMessage> messenger) {
        final Traversal.Admin<Vertex, Edge> edgeTraversal = this.edgeTraversal.getPure();
        edgeTraversal.addStart(edgeTraversal.getTraverserGenerator().generate(vertex, edgeTraversal.getStartStep(), 1));

        double tmin = Double.MAX_VALUE;
        Edge edgeMin = null;

        while (edgeTraversal.hasNext()) {
            final Edge edge = edgeTraversal.next();

            final AABB neighborAABB = (AABB) edge.property("aabb").value();

            Point3d p = new Point3d();

            if (neighborAABB.intersectRay(photon.ray(), p)) {
                if ((Math.abs(photon.ray().getTNear() - photon.exitPointT()) < EPSILON)
                        && (photon.ray().getTNear() < photon.ray().getTFar())
                        && (photon.ray().getTNear() < tmin)) {
                    tmin = photon.ray().getTNear();
                    edgeMin = edge;
                }
            }
        }

        if (edgeMin != null) {
            Vertex otherV = edgeMin.inVertex();

            if (otherV.equals(vertex))
                otherV = edgeMin.outVertex();

            messenger.sendMessage(MessageScope.Global.of(otherV), photon);

            return true;
        }

        return false;
    }

    private void collectPhoton(PhotonMessage photon, final Memory memory) {
        final List<Pair<Pair<Integer, Integer>, Color3f>> result;

        result = new ArrayList<>();

        result.add(new Pair<>(
                new Pair<>(photon.rayNumber() % width, photon.rayNumber() / width),
                new Color3f(photon.color())
        ));

        memory.add(IMAGE_COLORS_OUT_PARAM_STRING, result);
    }

    static class PhotonMessage {

        private String id;
        private Ray ray;
        private Color3f accColor;
        private Double accOpacity;
        private Integer rayNumber;
        Double lastSamplePointT;
        Double exitPointT;

        private final double EPSILON = 1e-6;

        private PhotonMessage(String id, Ray ray, Integer number) {
            this.id = id;
            this.rayNumber = number;
            this.accOpacity = 0.0;
            this.accColor = new Color3f(0.0f, 0.0f, 0.0f);
            this.lastSamplePointT = null;
            this.exitPointT = null;
            this.ray = new Ray(ray.getP(), ray.getD());
        }

        static PhotonMessage of(String id, Ray ray, Integer number) {
            return new PhotonMessage(id, ray, number);
        }

        Integer rayNumber() {
            return rayNumber;
        }

        Color3f color() {
            return accColor;
        }

        String id() {
            return id;
        }

        Ray ray() {
            return ray;
        }

        Double lastSamplePointT() {
            return lastSamplePointT;
        }

        Double exitPointT() {
            return exitPointT;
        }

        void setAccColor(Color3f color) {
            accColor.set(color);
        }

        void setLastSamplePointT(Double t) {
            if (lastSamplePointT == null)
                lastSamplePointT = new Double(t);
            else
                lastSamplePointT = t;
        }

        void setExitPointT(Double t) {
            if (exitPointT == null)
                exitPointT = new Double(t);
            else
                exitPointT = t;
        }

        void accumulateCT(Color3f color, Double opacity) {
            accColor.scaleAdd(1.0f - accOpacity.floatValue(), color, accColor);
            accColor.clampMax(1.0f);

            accOpacity = Math.min(accOpacity + (1.0 - accOpacity) * opacity, 1.0);
        }

        boolean isInsignificant() {
            return (Math.abs(1.0 - accOpacity) <= EPSILON) | accColor.epsilonEquals(new Point3f(1.0f, 1.0f, 1.0f), (float) EPSILON);
        }

        @Override
        public String toString() {
            return "Photon msg " + id() + ": ray = " + ray + ", pixel coordinates = " + rayNumber + ", accumulated color"
                    + accColor + ", accumulated opacity " + accOpacity;
        }
    }

    //////////////////////////////

    public static Builder build() {
        return new Builder();
    }

    @SuppressWarnings("WeakerAccess")
    public static final class Builder extends AbstractVertexProgramBuilder<Builder> {

        private Builder() {
            super(DVRVertexProgram.class);
        }

        public Builder height(final Integer height) {
            if (null != height)
                this.configuration.setProperty(HEIGHT_IN_PARAM_STRING, height);
            else
                this.configuration.clearProperty(HEIGHT_IN_PARAM_STRING);

            return this;
        }

        public Builder width(final Integer width) {
            if (null != width)
                this.configuration.setProperty(WIDTH_IN_PARAM_STRING, width);
            else
                this.configuration.clearProperty(WIDTH_IN_PARAM_STRING);

            return this;
        }

        public Builder P1(final Point3d P1) {
            if (null != P1)
                this.configuration.setProperty(P1_IN_PARAM_STRING, P1);
            else
                this.configuration.clearProperty(P1_IN_PARAM_STRING);

            return this;
        }

        public Builder origin(final Point3d origin) {
            if (null != origin)
                this.configuration.setProperty(ORIGIN_IN_PARAM_STRING, origin);
            else
                this.configuration.clearProperty(ORIGIN_IN_PARAM_STRING);

            return this;
        }

        public Builder DX(final Vector3d DX) {
            if (null != DX)
                this.configuration.setProperty(DX_IN_PARAM_STRING, DX);
            else
                this.configuration.clearProperty(DX_IN_PARAM_STRING);

            return this;
        }

        public Builder DY(final Vector3d DY) {
            if (null != DY)
                this.configuration.setProperty(DY_IN_PARAM_STRING, DY);
            else
                this.configuration.clearProperty(DY_IN_PARAM_STRING);

            return this;
        }

        public Builder samplingStep(final Number step) {
            if (null != step)
                this.configuration.setProperty(SAMPLING_STEP_IN_PARAM_STRING, step);
            else
                this.configuration.clearProperty(SAMPLING_STEP_IN_PARAM_STRING);
            return this;
        }

        public Builder pixelDistance(final Number distance) {
            if (null != distance)
                this.configuration.setProperty(PIXEL_DISTANCE_IN_PARAM_STRING, distance);
            else
                this.configuration.clearProperty(PIXEL_DISTANCE_IN_PARAM_STRING);
            return this;
        }

        public Builder sliceDistance(final Number distance) {
            if (null != distance)
                this.configuration.setProperty(SLICE_DISTANCE_IN_PARAM_STRING, distance);
            else
                this.configuration.clearProperty(SLICE_DISTANCE_IN_PARAM_STRING);
            return this;
        }

        public Builder sceneBox(final AABB sceneBox) {
            if (null != sceneBox)
                this.configuration.setProperty(SCENE_BOX_IN_PARAM_STRING, sceneBox);
            else
                this.configuration.clearProperty(SCENE_BOX_IN_PARAM_STRING);

            return this;
        }

    }

    ////////////////////////////

    @Override
    public Features getFeatures() {
        return new Features() {
            @Override
            public boolean requiresGlobalMessageScopes() {
                return true;
            }

            @Override
            public boolean requiresVertexPropertyAddition() {
                return true;
            }
        };
    }
}