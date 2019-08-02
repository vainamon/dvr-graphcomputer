package ru.sfedu.dvr;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import javax.vecmath.Point4i;
import javax.vecmath.Point3d;

import java.io.Serializable;

public class EquidistantVolumeCell implements Serializable {
    private static final Logger LOGGER = LoggerFactory.getLogger(EquidistantVolumeCell.class);

    private Point3d origin;
    private double pixelDistance = 1.0;
    private double sliceDistance = 1.0;
    // V000 -> V100 -> V010 -> V110
    private Point4i fIs;
    // V001 -> V101 -> V011 -> V111
    private Point4i bIs;

    public EquidistantVolumeCell(Point3d origin, Point4i frontIntensities, Point4i backIntensities) {
        this.origin = origin;
        this.fIs = frontIntensities;
        this.bIs = backIntensities;
    }

    public AABB getAABB() {
        return new AABB(origin, new Point3d(origin.getX() + pixelDistance,
                origin.getY() + pixelDistance, origin.getZ() + sliceDistance));
    }

    public double getTrilinearInterpolation(Point3d p) {
        double xWeight = Math.abs((p.getX() - origin.getX()) / pixelDistance);
        double yWeight = Math.abs((p.getY() - origin.getY()) / pixelDistance);
        double zWeight = Math.abs((p.getZ() - origin.getZ()) / sliceDistance);

        return  ((fIs.getX() * (1 - xWeight) + fIs.getY() * xWeight) * (1 - yWeight)
                + (fIs.getZ() * (1 - xWeight) + fIs.getW() * xWeight) * yWeight) * (1 - zWeight)
                + ((bIs.getX() * (1 - xWeight) + bIs.getY() * xWeight) * (1 - yWeight)
                + (bIs.getZ() * (1 - xWeight) + bIs.getW() * xWeight) * yWeight) * zWeight;
    }

    @Override
    public boolean equals(Object o) {
        if (o == this) {
            return true;
        }

        if (!(o instanceof EquidistantVolumeCell)) {
            return false;
        }

        EquidistantVolumeCell c = (EquidistantVolumeCell) o;

        return origin.equals(c.origin);
    }

    @Override
    public String toString() {
        return "EquCell: origin = " + origin;
    }
}
