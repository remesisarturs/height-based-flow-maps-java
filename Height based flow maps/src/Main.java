public class Main {


    public static void main(String[] args) {

        System.out.println("test");

        //Cell cell = new Cell();

        int nr_of_rows = 100;
        int nr_of_columns = 100;

        Cell[][] grid = initialize_grid(nr_of_rows, nr_of_columns);
        


    }

    public static Cell[][] initialize_grid(int number_of_rows, int number_of_columns){

        Cell[][] grid = new Cell[number_of_rows][number_of_columns];

        for (int i = 0; i < number_of_rows; i++) {
            for (int j = 0; j < number_of_columns; j++) {

                Cell cell = new Cell();

                cell.x = 0;
                cell.y = 0;
                cell.flow_direction = 0;
                cell.height = 0;

                cell.temp = Integer.toString(i) + " " + Integer.toString(j);

                grid[i][j] = cell;

                //System.out.println(grid[i][j].temp);
                //grid[i][j]
            }
        }

        return grid;
    }

}
