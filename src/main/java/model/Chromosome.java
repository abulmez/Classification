package model;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;


public class Chromosome implements Serializable {
    private ArrayList<Double> values;

    public Chromosome(Double... args){
        values = new ArrayList<>();
        values.addAll(Arrays.asList(args));
    }

    public ArrayList<Double> getValues() {
        return values;
    }

    public Double getValue(Integer position){
        return values.get(position);
    }

    public void setValue(Integer position, Double value) {
        values.set(position,value);
    }
}
