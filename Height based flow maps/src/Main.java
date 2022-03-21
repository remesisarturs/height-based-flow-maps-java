import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.util.*;

public class Main {


    public static void main(String[] args) throws FileNotFoundException {

        System.out.println("test");

        //Cell cell = new Cell();

        int nr_of_rows = 100;
        int nr_of_columns = 100;

        Cell[][] grid = initialize_grid(nr_of_rows, nr_of_columns);

        int source_x = 25;
        int source_y = 25;

        ArrayList items = read_input();

        Iterator it = items.iterator();

        while (it.hasNext()) {

            System.out.println(it.next());

        }


    }

    public static ArrayList read_input() throws FileNotFoundException {

        //parsing a CSV file into Scanner class constructor
        //System.out.println("Working Directory = " + System.getProperty("user.dir"));

        Scanner sc = new Scanner(new File("./input/1_s_2_t.csv"));
        sc.useDelimiter(",");   //sets the delimiter pattern

        ArrayList items = new ArrayList();

        while (sc.hasNext()) {
            //System.out.print(sc.next());

            items.add(sc.next());

        }
        sc.close();  //closes the scanner

        return items;
    }

    public static Cell[][] initialize_grid(int number_of_rows, int number_of_columns) {

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
