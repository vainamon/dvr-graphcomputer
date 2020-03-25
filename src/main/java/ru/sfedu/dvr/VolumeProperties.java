package ru.sfedu.dvr;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.Serializable;

public class VolumeProperties implements Serializable {
    private static final Logger LOGGER = LoggerFactory.getLogger(VolumeProperties.class);

    public Boolean isEmpty;
    public Double normalizedDeviation;

    public VolumeProperties() {
        this.isEmpty = false;
        this.normalizedDeviation = 0d;
    }

    @Override
    public String toString() {
        return "VolumeProperties: isEmpty = " + isEmpty + "; deviation = " + normalizedDeviation;
    }
}
