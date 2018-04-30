package model;

import java.util.ArrayList;

public class Matrix <E> {

    ArrayList<ArrayList<E>> matrix;
    Integer rows,columns;

    public Matrix(Integer rows, Integer columns, E initElem){
        this.columns = columns;
        this.rows = rows;
        matrix = new ArrayList<>();
        initMatrix(initElem);
    }

    private void initMatrix(E elem){
        for(int i=0;i<rows;i++){
            ArrayList<E> row = new ArrayList<>();
            for(int j=0;j<columns;j++){
                row.add(elem);
            }
            matrix.add(row);
        }
    }

    public E get(Integer... args){
        if(args.length == 1)
            return matrix.get(0).get(args[0]);
        return matrix.get(args[0]).get(args[1]);
    }


    public void set(E value,Integer... args){
        if(args.length == 1)
            matrix.get(0).set(args[0],value);
        else
            matrix.get(args[0]).set(args[1],value);
    }

    public Integer getRows() {
        return rows;
    }

    public Integer getColumns() {
        return columns;
    }

    public void setRows(Integer rows) {
        this.rows = rows;
    }
}
