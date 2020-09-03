package ru.sfedu.dvr;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import javax.vecmath.Point3d;
import javax.vecmath.Vector3d;

import java.io.Serializable;

public class AABB implements Serializable {
    private static final Logger LOGGER = LoggerFactory.getLogger(AABB.class);

    private final double EPSILON = 1e-6;

    private Point3d min;
    private Point3d max;

    private AABB()
    {
        // for kryo
    }

    public AABB(Point3d min, Point3d max) {
        this.min = min;
        this.max = max;
    }

    public Point3d min() {
        return this.min;
    }

    public Point3d max() {
        return this.max;
    }

    public boolean intersectRay(Ray r, Point3d q) {
        double tmin = 0.0d;
        double tmax = Double.MAX_VALUE;

        double d[] = new double[3];
        r.getD().get(d);

        double p[] = new double[3];
        r.getP().get(p);

        double aabbMin[] = new double[3];
        min.get(aabbMin);

        double aabbMax[] = new double[3];
        max.get(aabbMax);

        for (int i = 0; i < 3; i++) {
            if (Math.abs(d[i]) < EPSILON) {
// Ray is parallel to slab. No hit if origin not within slab
                if (p[i] < aabbMin[i] || p[i] > aabbMax[i]) return false;
            } else {
// Compute intersection t value of ray with near and far plane of slab
                double ood = 1.0d / d[i];
                double t1 = (aabbMin[i] - p[i]) * ood;
                double t2 = (aabbMax[i] - p[i]) * ood;

                double temp;
// Make t1 be intersection with near plane, t2 with far plane
                if (t1 > t2) {
                    temp = t1;
                    t1 = t2;
                    t2 = temp;
                }
// Compute the intersection of slab intersection intervals
                if (t1 > tmin) tmin = t1;
                if (t2 < tmax) tmax = t2;
// Exit with no collision as soon as slab intersection becomes empty
                //if ((tmin > tmax) || (Math.abs(tmin - tmax) < EPSILON) || (tmax < 0)) return false;
                if (tmin > tmax) return false;
            }
        }

// Ray intersects all 3 slabs. Return point (q) and intersection t value (tmin)
        q.scaleAdd(tmin, new Vector3d(r.getD()), new Point3d(r.getP()));

        r.setTNear(tmin);
        r.setTFar(tmax);

        return true;
    }

    public double sqDistanceToPoint(Point3d q) {
        double sqDist = 0.0f;

        double v[] = new double[3];
        q.get(v);

        double aabbMin[] = new double[3];
        min.get(aabbMin);

        double aabbMax[] = new double[3];
        max.get(aabbMax);

// For each axis count any excess distance outside box extents
        for (int i = 0; i < 3; i++) {
            if (v[i] < aabbMin[i]) sqDist += (aabbMin[i] - v[i]) * (aabbMin[i] - v[i]);
            if (v[i] > aabbMax[i]) sqDist += (v[i] - aabbMax[i]) * (v[i] - aabbMax[i]);
        }

        return sqDist;
    }

    @Override
    public String toString() {
        return "AABB: min = " + min + "; max = " + max;
    }
}
