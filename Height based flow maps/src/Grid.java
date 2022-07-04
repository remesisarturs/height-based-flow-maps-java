import java.util.ArrayList;
import java.util.Iterator;

public class Grid {

    public static Cell[][] initializeGrid(int nrOfRows, int nrOfColumns) {

        Cell[][] grid = new Cell[nrOfRows][nrOfColumns];

        for (int i = 0; i < nrOfRows; i++) {
            for (int j = 0; j < nrOfColumns; j++) {

                Cell cell = new Cell();

                cell.cellCol = i;
                cell.cellRow = j;
                cell.flowDirection = 0;
                cell.height = 0;

                cell.title = "";
                grid[i][j] = cell;
                cell.isObstacle = false;

                // i = columns
                // j = rows
                if (i == nrOfRows - 1 || j == nrOfColumns - 1 || i == 0 || j == 0) {
                    cell.title = "edge";
                }

                if (i == nrOfRows - 1) {
                    cell.title = "rightEdge";
                }

            }
        }

        return grid;
    }
    public static Bounds obtainBounds(ArrayList<Point> inputPoints) {

        Iterator it = inputPoints.iterator();

        double minX = inputPoints.get(0).col;
        double minY = inputPoints.get(0).row;
        double maxX = inputPoints.get(0).col;
        double maxY = inputPoints.get(0).row;

        while (it.hasNext()) {

            Point nextPoint = (Point) it.next();

            if (nextPoint.col < minX) {
                minX = nextPoint.col;
            }
            if (nextPoint.col > maxX) {
                maxX = nextPoint.col;
            }
            if (nextPoint.row < minY) {
                minY = nextPoint.row;
            }
            if (nextPoint.row > maxY) {
                maxY = nextPoint.row;
            }

        }

        Bounds result = new Bounds();
        result.maxX = maxX;
        result.minX = minX;
        result.minY = minY;
        result.maxY = maxY;

        return result;
    }

    public static void initializePointsInGrid(Cell[][] grid, ArrayList<Point> pointsList, String TARGET_NAME, int sourceX, int sourceY, Cell sourceCell) {

        Iterator it = pointsList.iterator();

        // setting the names of source/targets and setting height of source
        while (it.hasNext()) {

            Point point = (Point) it.next();

            grid[point.gridCol][point.gridRow].title = point.name;

            if (point.name.equals(TARGET_NAME)) {
                grid[point.gridCol][point.gridRow].height = -10000;

                sourceX = point.gridCol;
                sourceY = point.gridRow;

            }
        }

        sourceCell = grid[sourceX][sourceY];
    }

}
