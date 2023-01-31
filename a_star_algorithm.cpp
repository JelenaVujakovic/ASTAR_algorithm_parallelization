// A C++ Program to implement A* Search Algorithm
#include <bits/stdc++.h>
#include <mpi.h>
using namespace std;

#define ROW 10
#define COL 10

// Creating a shortcut for int, int pair type
typedef pair<int, int> Pair;
int grid[ROW][COL];
int subgrid1[ROW][COL];

// Create a closed list and initialise it to false which
// means that no cell has been included yet This closed
// list is implemented as a boolean 2D array
bool closedList[ROW][COL];
// Creating a shortcut for pair<int, pair<int, int>> type
typedef pair<double, pair<int, int> > pPair;
int ret = 0;
// We set this boolean value as false as initially
// the destination is not reached.
bool foundDest = false;

// A structure to hold the necessary parameters
struct cell {
	// Row and Column index of its parent
	// Note that 0 <= i <= ROW-1 & 0 <= j <= COL-1
	int parent_i, parent_j;
	// f = g + h
	double f, g;
};
// Declare a 2D array of structure to hold the details
// of that cell
    cell cellDetails[ROW][COL];

    /*
     Create an open list having information as-
     <f, <i, j>>
     where f = g + h,
     and i, j are the row and column index of that cell
     Note that 0 <= i <= ROW-1 & 0 <= j <= COL-1
     This open list is implemented as a set of pair of
     pair.*/
    set<pPair> openList;

    int* foundFixedPoint = new int;
 
int returnSize(char* fname){
    FILE* f = fopen(fname , "r");
    int dim = 0;
    char tmp;
    while (fscanf(f, "%c", &tmp) >= EOF && tmp >= '0'){
        dim++;
    }
    fclose(f);

    return dim;
}

int* loadMat(const char* fname, int dim)
{
    FILE* f = fopen(fname, "r");
    char* it = new char[dim*dim];
    static int src_dest[4]={0};
	int i=0;
	int j=0;
    while (fscanf(f, "%c", it) != EOF) {
        if(*it != '\n' && *it != '\r') 
		{
			if(*it == 'S'){
				src_dest[0] = i;
				src_dest[1] = j;
			}
			if(*it == 'E'){
				src_dest[2] = i;
				src_dest[3] = j;
			}

            if(*it == 'S' || *it=='E'){
				grid[i][j]='0';
			}else{
				grid[i][j]= *it;
			}
			it++;
			j++;
			if(j == COL){
				j=0;
				i++;
			}
		}
    }
    fclose(f);
    return src_dest;
}

// A Utility Function to check whether given cell (row, col)
// is a valid cell or not.
bool isValid(int row, int col)
{
	// Returns true if row number and column number
	// is in range
	return (row >= 0) && (row < ROW) && (col >= 0)
		&& (col < COL);
}

// A Utility Function to check whether the given cell is
// blocke``d or not
bool isUnBlocked(int grid[][COL], int row, int col)
{
	// Returns true if the cell is not blocked else false
	if (grid[row][col] == 0  || grid[row][col] == 48){
		return (true);
	}
	else{
		return (false);
	}
}

// A Utility Function to check whether destination cell has
// been reached or not
bool isDestination(int row, int col, Pair dest)
{
	if (row == dest.first && col == dest.second)
		return (true);
	else
		return (false);
}


// A Utility Function to trace the path from the source
// to destination
// A Utility Function to trace the path from the source
// to destination
void tracePath(cell cellDetails[][COL], Pair dest, Pair src, int prank, int is_in_range,int csize)
{  
        int row = dest.first;
        int col = dest.second;
        int i = 1;

        stack<Pair> Path;
        FILE *fptr;

        if ((fptr = fopen("res.txt","wb+")) == NULL){
            printf("Error! opening file");
            // Program exits if the file pointer returns NULL.
            exit(1);
        }
    
        while (!(cellDetails[row][col].parent_i == row && cellDetails[row][col].parent_j == col)) {
            Path.push(make_pair(row, col));
            int temp_row = cellDetails[row][col].parent_i;
            int temp_col = cellDetails[row][col].parent_j;
            row = temp_row;
            col = temp_col;
            
        }
        
        Path.push(make_pair(row, col));
        while (!Path.empty()) {
            pair<int, int> p = Path.top();  
            Path.pop();
            if((prank == csize-1) && (is_in_range == 0) || (prank == is_in_range)){
                if(i==1){
                    if(is_in_range==prank){
                    printf("Process %d found path\n",prank);
                    }
                    printf("The Path is ");
                    i=0;
                }
                printf("-> (%d,%d) ", p.first, p.second);
                grid[p.first][p.second]='x';
            }
        }
        grid[src.first][src.second]='S';
        grid[dest.first][dest.second]='E';
        for(int i=0;i<ROW;i++){
            for(int j=0;j<COL;j++){
                if((grid[i][j] == 48 || grid[i][j] == 49 || grid[i][j] == 0 || grid[i][j] == 1) && grid[i][j] != 'x'){
                    grid[i][j] = ' ';
                }
                fputc(grid[i][j], fptr);
                if(j==COL-1)
                    fputc('\n',fptr);
            }
        }
        fclose(fptr);
        return;
}
void log_result(int prank, int is_in_range, Pair src,Pair dest, int csize){
    if(prank == 1){
        *foundFixedPoint = 1;
        if(is_in_range != 1){
            MPI_Send(foundFixedPoint,1,MPI_INT,2,0,MPI_COMM_WORLD);
        }
    }else if(prank > 1){
        if(is_in_range != prank){
            MPI_Recv(foundFixedPoint,1,MPI_INT,prank-1,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            if(prank < csize-1){
                MPI_Send(foundFixedPoint,1,MPI_INT,prank+1,0,MPI_COMM_WORLD);
            }
        }
    }
}  
// Only process this cell if this is a valid one
bool traverse_matrix(int k, int i, int r, int j,Pair src, Pair dest, int prank, int csize, int is_in_range,double gNew,double fNew){
    if (isValid(k, r) == true) {
        // If the destination cell is the same as the
        // current successor
        if (isDestination(k, r, dest) == true) {
            // Set the Parent of the destination cell
            cellDetails[k][r].parent_i = i;
            cellDetails[k][r].parent_j = j;
            log_result(prank,is_in_range,src,dest, csize);
            tracePath(cellDetails, dest, src,prank,is_in_range,csize);
            foundDest = true;
            return true;
        }
        // If the successor is already on the closed
        // list or if it is blocked, then ignore it.
        // Else do the following
        else if (closedList[k][r] == false
                 && isUnBlocked(grid, k, r)
                        == true) {
            gNew = cellDetails[i][j].g + 1.0;
            //hNew = calculateHValue(k, j, dest);
            fNew = gNew;
            // If it isn’t on the open list, add it to
            // the open list. Make the current square
            // the parent of this square. Record the
            // f, g, and h costs of the square cell
            //                OR
            // If it is on the open list already, check
            // to see if this path to that square is
            // better, using 'f' cost as the measure.
            if (cellDetails[k][r].f == FLT_MAX
                || cellDetails[k][r].f > fNew) {
                openList.insert(make_pair(
                    fNew, make_pair(k, r)));    
                // Update the details of this cell
                cellDetails[k][r].f = fNew;
                cellDetails[k][r].g = gNew;
                cellDetails[k][r].parent_i = i;
                cellDetails[k][r].parent_j = j;
            }
        }
    }
    return false;
}

// A Function to find the shortest path between
// a given source cell to a destination cell according
// to A* Search Algorithm
void aStarSearch(int grid[][COL], Pair src, Pair dest, int prank, int csize, int is_in_range)
{
    int i, j;
 
    // Initialising the parameters of the starting node
    i = src.first, j = src.second;
    cellDetails[i][j].f = 0.0;
    cellDetails[i][j].g = 0.0;
    cellDetails[i][j].parent_i = i;
    cellDetails[i][j].parent_j = j;

     // Put the starting cell on the open list and set its
    // 'f' as 0
    openList.insert(make_pair(0.0, make_pair(i, j)));
 
 
    while (!openList.empty()) {
        pPair p = *openList.begin();
 
        // Remove this vertex from the open list
        openList.erase(openList.begin());
 
        // Add this vertex to the closed list
        i = p.second.first;
        j = p.second.second;
        closedList[i][j] = true;
        /*
         Generating all the 8 successor of this cell
 
             N.W   N   N.E
               \   |   /
                \  |  /
             W----Cell----E
                  / | \
                /   |  \
             S.W    S   S.E
 
         Cell-->Popped Cell (i, j)
         N -->  North       (i-1, j)
         S -->  South       (i+1, j)
         E -->  East        (i, j+1)
         W -->  West           (i, j-1)
         N.E--> North-East  (i-1, j+1)
         N.W--> North-West  (i-1, j-1)
         S.E--> South-East  (i+1, j+1)
         S.W--> South-West  (i+1, j-1)*/
 
        // To store the 'g', 'h' and 'f' of the 8 successors
         // Only process this cell if this is a valid one
        double gNew, fNew;
        //----------- 1st Successor (North) ------------
        if(true == traverse_matrix(i-1,i,j,j,src, dest, prank, csize, is_in_range,gNew,fNew)){
            break;
        }
        //----------- 2nd Successor (South) ------------
        if(true == traverse_matrix(i+1,i,j,j,src, dest, prank, csize, is_in_range,gNew,fNew)){
            break;
        }
        //----------- 3rd Successor (East) ------------
        if(true == traverse_matrix(i,i,j+1,j,src, dest, prank, csize, is_in_range,gNew,fNew)){
            break;
        }
        //----------- 4th Successor (West) ------------
        if(true == traverse_matrix(i,i,j-1,j,src, dest, prank, csize, is_in_range,gNew,fNew)){
            break;
        } 
    }
    // When the destination cell is not found and the open
    // list is empty, then we conclude that we failed to
    // reach the destination cell. This may happen when the
    // there is no way to destination cell (due to
    // blockages)
    if (foundDest == false)
        printf("Failed to find the Destination Cell for process %d \n",prank);
 
    return;
}

int check_validity(Pair src, Pair dest, int prank){
            // If the source is out of range
            if (isValid(src.first, src.second) == false) {
                printf("Source is invalid\n");
                return 0;
            }
        
            // If the destination is out of range
            if (isValid(dest.first, dest.second) == false) {
                printf("Destination is invalid\n");
                return 0;
            }
        
            // Either the source or the destination is blocked
            if (isUnBlocked(grid, src.first, src.second) == false) {
                printf("Source is blocked\n");
                return 0;
            }
            if (isUnBlocked(grid, dest.first, dest.second) == false) {
                printf("Destination is blocked\n");
                return 0;
            }
            
            // If the destination cell is the same as source cell
            if (isDestination(src.first, src.second, dest)
                == true) {
                printf("We are already at the destination %d\n",prank);
                return 0;
            }

            return 1;
}
int check_matrix(int csize, int dim, int prank,int src,int dst)
{
    int dim_submatrix;

    dim_submatrix = (int)(dim/(csize-1));

    //Check if S and E are in submatrix of 1st process
    if(src < dim_submatrix && dst < dim_submatrix)
    {
        return 1;
    }

    if(csize >= 3)
    {
        //Check if S and E are in submatrix of 2nd process
        if(src >= dim_submatrix && dst >= dim_submatrix){
            return 2;
        }
    }
    if(csize >= 4)
    {
        //Check if S and E are in submatrix of 2nd process
        if(((src >= dim_submatrix && src < dim_submatrix+1) && (dst >= dim_submatrix && dst < dim_submatrix+1))){
            return 2;
            //Check if S and E are in submatrix of 3rd process
        }else if((src >= dim_submatrix+1 && dst >= dim_submatrix+1)){
            return 3;
        }
    }

    if(csize >= 5)
    {
        //Check if S and E are in submatrix of 2nd process
        if(((src >= dim_submatrix && src < dim_submatrix+2) && (dst >= dim_submatrix && dst < dim_submatrix+2))){
            return 2;
        //Check if S and E are in submatrix of 3rd process
        }else if(((src >= dim_submatrix+2 && src < dim_submatrix+4) && (dst >= dim_submatrix+2 && dst < dim_submatrix+4))){
            return 3;
        }else if((src >= dim_submatrix+4 && dst >= dim_submatrix+2)){
            return 4;
        }
    }
    if(csize >= 6)
    {
        //Check if S and E are in submatrix of 2nd process
        if(((src >= dim_submatrix && src < dim_submatrix+2) && (dst >= dim_submatrix && dst < dim_submatrix+2))){
            return 2;
        //Check if S and E are in submatrix of 3rd process
        }else if(((src >= dim_submatrix+2 && src < dim_submatrix+4 ) && (dst >= dim_submatrix+2 && dst < dim_submatrix+4))){
            return 3;
        }else if(((src >= dim_submatrix+4 && src < dim_submatrix+6) && (dst >= 6 && dst < dim_submatrix+6))){
            return 4;
        }
        else if(src >= dim_submatrix+8 && (dst >= dim_submatrix+8)){
            return 5;
        }
    }
    return 0;
}

void start_a_star_on_process(int prank, int src_x, int src_y,int dst_x,int dst_y,int is_in_range, int dim, int csize,int fixed_src_x,int fixed_src_y,int fixed_dst_x,int fixed_dst_y){
    Pair src;
    Pair dest_fixed_point;

    if(src_x <= fixed_src_x)
    {
        src = make_pair(src_x,src_y);
    }
    else{
        src = make_pair(fixed_src_x,fixed_src_y);
    }
    // Destination is the left-most top-most corner
    if(dst_x <= fixed_src_x){
        // Destination is the left-most top-most corner
        dest_fixed_point = make_pair(dst_x, dst_y);
    }else{
        dest_fixed_point = make_pair(fixed_dst_x,fixed_dst_y);
    }
    check_validity(src,dest_fixed_point, prank);
    for(int i = 0;i<COL/2;i++){
        for(int j=0;j<COL;j++){
            subgrid1[i][j] = grid[i][j]; 
        }
    }
    if((is_in_range == prank) || (is_in_range == 0)){
        aStarSearch(subgrid1, src, dest_fixed_point,prank,csize,is_in_range);
    }
}

// Driver program to test above function
int main(int argc, char* argv[])
{
    int csize;
    int prank;
    int* bcast_array;
    int dim;
    FILE *fptr;

    if ((fptr = fopen("performance_metrics.txt","a+")) == NULL){
        printf("Error! opening file");
        // Program exits if the file pointer returns NULL.
        exit(1);
    }
    
    MPI_Init (NULL,NULL);
    MPI_Comm_size (MPI_COMM_WORLD,&csize);
    MPI_Comm_rank (MPI_COMM_WORLD,&prank);

    if(prank == 0)
    {
        dim = returnSize(argv[1]);
        //Load matrix
        bcast_array = loadMat(argv[1], dim);
        bcast_array[5] = dim;

        if(csize > 1)
            bcast_array[4] = check_matrix(csize,dim,prank,bcast_array[0],bcast_array[2]);
    }
    else{
        bcast_array = new int[7];
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(bcast_array,7,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    double s = MPI_Wtime();

    dim = bcast_array[5];
    Pair src;
    Pair dest;
    *foundFixedPoint = 0;
    memset(closedList, false, sizeof(closedList));
    for(int i = 0; i < ROW; i++){
        for (int j = 0; j < COL; j++) {
            cellDetails[i][j].f = FLT_MAX;
            cellDetails[i][j].g = FLT_MAX;
            cellDetails[i][j].parent_i = -1;
            cellDetails[i][j].parent_j = -1;
        }
    }
    //Sequential code
    if(csize < 3){
        if((prank==1 && csize == 2) || (csize == 1 && prank == 0)){
            start_a_star_on_process(prank,bcast_array[0],bcast_array[1],bcast_array[2],bcast_array[3],bcast_array[4],dim,csize,dim,dim,dim,dim);
        }
    }
    else if(csize == 3){
        if(prank == 1){
            start_a_star_on_process(prank,bcast_array[0],bcast_array[1],bcast_array[2],bcast_array[3],bcast_array[4],dim,csize,prank-1,prank-1,(dim/2)-1,dim-1);
        } 
        if(prank == 2){
            start_a_star_on_process(prank,bcast_array[0],bcast_array[1],bcast_array[2],bcast_array[3],bcast_array[4],dim,csize,(dim/2)-1,dim-1,dim-1,dim-1);
        }
    }else if(csize == 4){
        if(prank == 1){
            start_a_star_on_process(prank,bcast_array[0],bcast_array[1],bcast_array[2],bcast_array[3],bcast_array[4],dim,csize,prank-1,prank-1,dim/5,dim-1);
        }else if(prank == 2){
            start_a_star_on_process(prank,bcast_array[0],bcast_array[1],bcast_array[2],bcast_array[3],bcast_array[4],dim,csize,dim/5,dim-1,(dim/2)-1,dim-1);
        }else if(prank == 3){
            start_a_star_on_process(prank,bcast_array[0],bcast_array[1],bcast_array[2],bcast_array[3],bcast_array[4],dim,csize,(dim/2)-1,dim-1,dim-1,dim-1);
        }
    }else if(csize == 5){
        if(prank == 1){
            start_a_star_on_process(prank,bcast_array[0],bcast_array[1],bcast_array[2],bcast_array[3],bcast_array[4],dim,csize,prank-1,prank-1,dim/5,dim-1);
        }else if(prank == 2){
            start_a_star_on_process(prank,bcast_array[0],bcast_array[1],bcast_array[2],bcast_array[3],bcast_array[4],dim,csize,dim/5,dim-1,(dim/2)-1,dim-1);
        }else if(prank == 3){
            start_a_star_on_process(prank,bcast_array[0],bcast_array[1],bcast_array[2],bcast_array[3],bcast_array[4],dim,csize,(dim/2)-1,dim-1,(dim/2)+1,dim-1);
        }else if(prank == 4){
            start_a_star_on_process(prank,bcast_array[0],bcast_array[1],bcast_array[2],bcast_array[3],bcast_array[4],dim,csize,(dim/2)+1,dim-1,dim-1,dim-1);
        }
    }else if(csize == 6){
        if(prank == 1){
            start_a_star_on_process(prank,bcast_array[0],bcast_array[1],bcast_array[2],bcast_array[3],bcast_array[4],dim,csize,prank-1,prank-1,dim/5,dim-1);
        }else if(prank == 2){
            start_a_star_on_process(prank,bcast_array[0],bcast_array[1],bcast_array[2],bcast_array[3],bcast_array[4],dim,csize,dim/5,dim-1,(dim/2)-1,dim-1);
        }else if(prank == 3){
            start_a_star_on_process(prank,bcast_array[0],bcast_array[1],bcast_array[2],bcast_array[3],bcast_array[4],dim,csize,(dim/2)-1,dim-1,(dim/2)+1,dim-1);
        }else if(prank == 4){
            start_a_star_on_process(prank,bcast_array[0],bcast_array[1],bcast_array[2],bcast_array[3],bcast_array[4],dim,csize,(dim/2)+1,dim-1,dim-2,dim-1);
        }else if(prank == 5){
            start_a_star_on_process(prank,bcast_array[0],bcast_array[1],bcast_array[2],bcast_array[3],bcast_array[4],dim,csize,dim-2,dim-1,dim-1,dim-1);
        }
    }
    double e = MPI_Wtime();
    if(prank==0)
        printf("\nTime elapsed: %lf\n",e-s);
    MPI_Finalize ();
    if(prank==0){
        fprintf(fptr,"Time elapsed: %lf\t",e-s);
        fprintf(fptr,"Number of proceses = %d\t",csize);
        fprintf(fptr,"Input matrix size = %dx%d\n",ROW,COL);
    }
	return (0);
}

