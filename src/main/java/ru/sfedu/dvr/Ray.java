package ru.sfedu.dvr;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import javax.vecmath.Point3d;
import javax.vecmath.Vector3d;

import java.io.Serializable;

public class Ray implements Serializable {
    private static final Logger LOGGER = LoggerFactory.getLogger(Ray.class);

    private Point3d P;
    private Vector3d d;
    double tNear = 0.0;
    double tFar = 0.0;

    public Ray(Point3d P, Vector3d d) {
        this.P = P;
        this.d = d;
    }

    public Vector3d getD() {
        return d;
    }

    public Point3d getP() {
        return P;
    }

    public void setTNear(double tNear) {
        this.tNear = tNear;
    }

    public double getTNear() {
        return this.tNear;
    }

    public void setTFar(double tFar) {
        this.tFar = tFar;
    }

    public double getTFar() {
        return this.tFar;
    }

    public Point3d getPointByT(double _t) {
        Vector3d _d = new Vector3d(d);
        _d.scale(_t);
        _d.add(P);

        return new Point3d(_d);
    }

    @Override
    public String toString() {
        return "Ray (P + t*d): " + P + " + t * " + d + "; tn = " + tNear + ", tf = " + tFar;
    }
}
