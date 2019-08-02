package ru.sfedu.dvr;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import javax.vecmath.Point3d;
import javax.vecmath.Tuple3d;
import javax.vecmath.Vector3d;

import java.io.Serializable;

public class Camera implements Serializable {
    private static final Logger LOGGER = LoggerFactory.getLogger(Camera.class);

    private Point3d eye;

    private Vector3d U;
    private Vector3d V;
    private Vector3d N;

    public Camera() {
        eye = new Point3d(0, 0 , 0);

        N = new Vector3d(0, 0, 1);
        U = new Vector3d(0, 1, 0);
        V = new Vector3d(1, 0, 0);
    }

    public Vector3d getN() {
        return  N;
    }

    public Tuple3d transform(Tuple3d v) {
        Vector3d eyeV = new Vector3d(eye.getX(), eye.getY(), eye.getZ());

        Tuple3d result;

        if (v instanceof Vector3d)
            result = new Vector3d();
        else
            result = new Point3d();

        final double x,y,z;

        x = v.getX() * U.getX() + v.getY() * U.getY()
                + v.getZ() * U.getZ() + (-eyeV.dot(U));

        y = v.getX() * V.getX() + v.getY() * V.getY()
                + v.getZ() * V.getZ() + (-eyeV.dot(V));

        z = v.getX() * N.getX() + v.getY() * N.getY()
                + v.getZ() * N.getZ() + (-eyeV.dot(N));

        result.set(x,y,z);

        return result;
    }

    void set(Point3d eye, Point3d look, Vector3d up)
    {
        this.eye.set(eye);

        Vector3d diff = new Vector3d(eye);

        diff.sub(look);

        N.set(diff);

        Vector3d cross = new Vector3d();

        cross.cross(up, N);

        U.set(cross);

        N.normalize();
        U.normalize();

        cross.cross(N, U);

        V.set(cross);
    }

    void slide(double _du, double _dv, double _dn)
    {
        final double x = eye.getX() + _du * U.getX() + _dv * V.getX() + _dn * N.getX();
        final double y = eye.getY() + _du * U.getY() + _dv * V.getY() + _dn * N.getY();
        final double z = eye.getZ() + _du * U.getZ() + _dv * V.getZ() + _dn * N.getZ();

        eye.set(x,y,z);
    }

    void roll(double _angle)
    {
        double cs = Math.cos(Math.PI / 180 * _angle);
        double sn = Math.sin(Math.PI / 180 * _angle);

        Vector3d tmp = new Vector3d(U);

        U.set(cs * tmp.getX() + sn * V.getX(),cs * tmp.getY() + sn * V.getY(),cs * tmp.getZ() + sn * V.getZ());
        V.set(-sn * tmp.getX() + cs * V.getX(),-sn * tmp.getY() + cs * V.getY(),-sn * tmp.getZ() + cs * V.getZ());
    }

    void pitch(double _angle)
    {
        double cs = Math.cos(Math.PI / 180 * _angle);
        double sn = Math.sin(Math.PI / 180 * _angle);

        Vector3d tmp = new Vector3d(N);

        N.set(cs * tmp.getX() + sn * V.getX(),cs * tmp.getY() + sn * V.getY(),cs * tmp.getZ() + sn * V.getZ());
        V.set(-sn * tmp.getX() + cs * V.getX(),-sn * tmp.getY() + cs * V.getY(),-sn * tmp.getZ() + cs * V.getZ());
    }

    void yaw(double _angle)
    {
        double cs = Math.cos(Math.PI / 180 * _angle);
        double sn = Math.sin(Math.PI / 180 * _angle);

        Vector3d tmp = new Vector3d(N);

        N.set(cs * tmp.getX() + sn * U.getX(),cs * tmp.getY() + sn * U.getY(),cs * tmp.getZ() + sn * U.getZ());
        U.set(-sn * tmp.getX() + cs * U.getX(),-sn * tmp.getY() + cs * U.getY(),-sn * tmp.getZ() + cs * U.getZ());
    }
}
