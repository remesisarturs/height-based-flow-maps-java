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

        ArrayList points_list = process_input(items);


        System.out.println();


    }

    public static ArrayList process_input(ArrayList items) {

        Iterator it = items.iterator();

        ArrayList point_list = new ArrayList();

        while (it.hasNext()) {

            String input_line = it.next().toString();
            List<String> components = Arrays.asList(input_line.split(";"));

            Point point = new Point();
            point.name = components.get(0);
            point.x = Float.parseFloat(components.get(1));
            point.y = Float.parseFloat(components.get(2));

            point_list.add(point);

        }

        return point_list;
    }

    public static ArrayList obtain_bounds(ArrayList input_points) {

//        for item in input_points:
//
//        if input_points.index(item) == 0:
//        minx = item[1]
//        miny = item[2]
//        maxx = item[1]
//        maxy = item[2]
//
//        if item[1] < minx:
//        minx = item[1]
//        if item[1] > maxx:
//        maxx = item[1]
//        if item[2] < miny:
//        miny = item[2]
//        if item[2] > maxy:
//        maxy = item[2]
//
//        return minx, maxx, miny, maxy

        Iterator it = input_points.iterator();

        while (it.hasNext()) {


        }
        return null;
    }


    public static ArrayList read_input() throws FileNotFoundException {

        //parsing a CSV file into Scanner class constructor
        //System.out.println("Working Directory = " + System.getProperty("user.dir"));

        Scanner sc = new Scanner(new File("./input/1_s_2_t.csv"));
        sc.useDelimiter("\n");   //sets the delimiter pattern

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
