package model;

import java.util.ArrayList;

public class Population {
    ArrayList<Chromosome> population;

    public  Population(){
        population = new ArrayList<>();
    }

    public void addChromosome(Chromosome c){
        population.add(c);
    }

    public ArrayList<Chromosome> getChromosomesList(){
        return population;
    }

    public Chromosome getChromosome(Integer position){
        return population.get(position);
    }

}
