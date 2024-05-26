#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <omp.h>
#include <unistd.h>
#include <time.h>
#include <math.h>

#define NUM_PORTS 5
#define FULL_THRESHOLD 3
#define TABLE_SIZE 5
#define CYCLE_TIME 5

#define REQUEST_TAG 0
#define AVAILABILITY_TAG 1
#define REPORT_TAG 2
#define NEARBY_NODES_TAG 3
#define TERMINATION_TAG 4
#define COMM_TIME_TAG 5

// Array for each charging station
typedef struct {
    int year, month, day, hour, minute, second;
    int availability;
} TableEntry;

typedef struct {
    TableEntry data[TABLE_SIZE];
    int count;
} Table;

typedef struct {
    int year, month, day, hour, minute, second;
    double start_time; 
    int coords[2];
    int availability;
    int neighbours[4];
    int nbr_availabilities[4];
    int numMessagesExchanged;
} Report;

typedef struct {
    double total_comm_time;
    int num_communications;
} CommTime;

MPI_Datatype create_mpi_report() {
    const int nitems = 12; // number of elements in the structure
    int blocklengths[12] = {1, 1, 1, 1, 1, 1, 1, 2, 1, 4, 4, 1}; // number of elements in each block
    MPI_Datatype types[12] = {MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_DOUBLE, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT}; // types of elements
    MPI_Datatype mpi_report;
    MPI_Aint offsets[12]; 

    offsets[0] = offsetof(Report, year);
    offsets[1] = offsetof(Report, month);
    offsets[2] = offsetof(Report, day);
    offsets[3] = offsetof(Report, hour);
    offsets[4] = offsetof(Report, minute);
    offsets[5] = offsetof(Report, second);
    offsets[6] = offsetof(Report, start_time);
    offsets[7] = offsetof(Report, coords);
    offsets[8] = offsetof(Report, availability);
    offsets[9] = offsetof(Report, neighbours);
    offsets[10] = offsetof(Report, nbr_availabilities);
    offsets[11] = offsetof(Report, numMessagesExchanged);

    MPI_Type_create_struct(nitems, blocklengths, offsets, types, &mpi_report);
    MPI_Type_commit(&mpi_report);

    return mpi_report;
}

MPI_Datatype create_mpi_commtime() {
    const int nitems = 2; // number of elements in the structure
    int blocklengths[2] = {1, 1}; // number of elements in each block
    MPI_Datatype types[2] = {MPI_DOUBLE, MPI_INT}; // types of elements
    MPI_Datatype mpi_commtime_type;
    MPI_Aint offsets[2];

    offsets[0] = offsetof(CommTime, total_comm_time);
    offsets[1] = offsetof(CommTime, num_communications);

    MPI_Type_create_struct(nitems, blocklengths, offsets, types, &mpi_commtime_type);
    MPI_Type_commit(&mpi_commtime_type);

    return mpi_commtime_type;
}


// Insert new entry at top and shift the rest down
void insertTable(Table *table, TableEntry entry) {
    for (int i = TABLE_SIZE - 1; i > 0; i--) {
        table->data[i] = table->data[i - 1];
    }
    table->data[0] = entry;
    if (table->count < TABLE_SIZE) {
        table->count++;
    }
}

// Remove the last entry from the table
TableEntry removeTable(Table *table) {
    TableEntry entry;
    if (table->count == 0) {
        entry.year = -1; // Represents an empty entry
        return entry;
    }
    entry = table->data[table->count - 1];
    table->count--;
    return entry;
}

int checkNodeAvailability(int ports[NUM_PORTS]) {
    int result = 0;
    for (int i = 0; i < NUM_PORTS; i++) {
        result += 1 - ports[i];
    }
    return result;
}

void write_report_to_file(const Report *report, const char *filename, int iteration, int node_rank) {
    // Open the file in append mode
    FILE *file = fopen(filename, "a"); 
    if (file == NULL) {
        perror("Error opening file");
        return;
    }

    // Write the report content to the file
    fprintf(file, "Report\n");
    fprintf(file, "Iteration: %d\n", iteration);
    fprintf(file, "Charging Node: %d\n", node_rank);
    fprintf(file, "Year: %d\n", report->year);
    fprintf(file, "Month: %d\n", report->month);
    fprintf(file, "Day: %d\n", report->day);
    fprintf(file, "Hour: %d\n", report->hour);
    fprintf(file, "Minute: %d\n", report->minute);
    fprintf(file, "Second: %d\n", report->second);
    fprintf(file, "Coords: (%d, %d)\n", report->coords[0], report->coords[1]);
    fprintf(file, "Availability: %d\n", report->availability);
    fprintf(file, "Neighbours: %d, %d, %d, %d\n", 
        report->neighbours[0], report->neighbours[1], report->neighbours[2], report->neighbours[3]);
    fprintf(file, "Neighbour Availabilities: %d, %d, %d, %d\n", 
        report->nbr_availabilities[0], report->nbr_availabilities[1], 
        report->nbr_availabilities[2], report->nbr_availabilities[3]);
    fprintf(file, "Num Messages Exchanged: %d\n", report->numMessagesExchanged);
    fprintf(file, "-----------------------------------------------------\n\n");

    // Close the file
    fclose(file);
}

int* check_neighbours_of_neighbours(Report node_reports[], int node_rank, int NUM_NODES, int iteration, int last_report[]) {
    // Get current time as epoch time
    time_t current_time;
    time(&current_time);

    // Allocate memory to store the IDs of nodes without recent reports
    int *nodes_without_reports = (int*)malloc(NUM_NODES * sizeof(int));

    // Retrieve the neighbours from the report
    int *neighbours = node_reports[node_rank].neighbours;

    // Loop through the neighbours
    for (int i = 0; i < 4; i++) { // Assuming maximum of 4 neighbors
        int neighbour_id = neighbours[i];
        if (neighbour_id >= 0 && neighbour_id < NUM_NODES) { // If the neighbour is valid
            // Access the neighbours of the neighbour
            int *neighbours_of_neighbour = node_reports[neighbour_id].neighbours;

            // Loop through the neighbours of the neighbour
            for (int j = 0; j < 4; j++) { // Assuming maximum of 4 neighbors
                int neighbour_of_neighbour_id = neighbours_of_neighbour[j];
                if (neighbour_of_neighbour_id >= 0 && neighbour_of_neighbour_id < NUM_NODES) { // Check if the neighbour of neighbour is valid
                    if (neighbour_of_neighbour_id != node_rank) {
                        if (iteration - last_report[neighbour_of_neighbour_id]) {
                        // Add this node's ID to the array
                        nodes_without_reports[neighbour_of_neighbour_id] = 1;
                    }
                    }

                }
            }
        }
    }

    return nodes_without_reports; // Return the array of nearby nodes
}


int runBaseStation(int NUM_NODES) {

    int terminate = 0; // terminates program

    Report received_report;
    Report node_reports[NUM_NODES];
    int last_report[NUM_NODES];
    MPI_Status received_report_status;
    int flag;
    int iteration = 1;

    double total_comm_time = 0;
    double report_end_time;
    double report_start_time;
    double comm_time;

    int num_reports = 0;

    // initialises values to -1
    for (int i = 0; i < NUM_NODES; ++i) {
        last_report[i] = -3;
    }

    MPI_Datatype mpi_report_type = create_mpi_report(); // creates mpi_report data type

    #pragma omp parallel num_threads(NUM_NODES + 1), private(received_report, received_report_status, flag, iteration, report_start_time, report_end_time, comm_time), shared(terminate, total_comm_time, num_reports)
    {
        int thread_id = omp_get_thread_num();

        if (thread_id >= 0 && thread_id <= NUM_NODES - 1) {

            iteration = 1;
            int node_rank = thread_id + 1;


            while (terminate == 0) {
                flag = 0; // flag to check if there is a message waiting
                #pragma omp critical
                {
                MPI_Iprobe(node_rank, REPORT_TAG, MPI_COMM_WORLD, &flag, MPI_STATUS_IGNORE);
                }
                if (flag) {
                    #pragma omp critical
                    {
                    MPI_Recv(&received_report, 1, mpi_report_type, node_rank, REPORT_TAG, MPI_COMM_WORLD, &received_report_status);

                    report_end_time = MPI_Wtime();

                    report_start_time = received_report.start_time; 

                    comm_time = report_end_time - report_start_time;
                    
                    total_comm_time += comm_time;

                    num_reports++;
                    }

                    printf("report start time: %f\n", report_start_time);
                    printf("report end time: %f\n", report_end_time);
                    printf("total comm time: %f\n", total_comm_time);
                    printf("num reports: %d\n", num_reports);

                    node_reports[thread_id] = received_report;
                    last_report[thread_id] = iteration;
                    int* nearby_nodes = check_neighbours_of_neighbours(node_reports, thread_id, NUM_NODES, iteration, last_report);
                    #pragma omp critical
                    {
                    write_report_to_file(&received_report, "log.txt", iteration, thread_id);
                    }

                    #pragma omp critical
                    {
                    MPI_Send(nearby_nodes, NUM_NODES, MPI_INT, node_rank, NEARBY_NODES_TAG, MPI_COMM_WORLD);
                    }
                }
            
            }
            iteration++;
            sleep(CYCLE_TIME);

        } else if (thread_id == NUM_NODES) {
            MPI_Datatype mpi_commtime_type = create_mpi_commtime(); // declare mpi_commtime_type

            printf("Press Enter to Terminate\n");
            int termination_message = 1;
            CommTime node_comm_time;
            double node_total_comm_time = 0;
            int node_num_comms = 0;
            double avg_node_comm_time;
            while (terminate == 0) {
                getchar();
                for (int i = 1; i < NUM_NODES + 1; ++i) { // Starting from 1 because 0 is the root itself
                    MPI_Send(&termination_message, 1, MPI_INT, i, TERMINATION_TAG, MPI_COMM_WORLD); // Send the termination message to each process
                    MPI_Recv(&node_comm_time, 1, mpi_commtime_type, i, COMM_TIME_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                    node_total_comm_time += node_comm_time.total_comm_time; // Add to the total communication time
                    node_num_comms += node_comm_time.num_communications;
                }

                // printf("total comm time: %f\n", total_comm_time);
                // printf("num reports: %d\n", num_reports);

                
                MPI_Type_free(&mpi_commtime_type);

                avg_node_comm_time = node_total_comm_time / node_num_comms;

                double avg_base_comm_time = total_comm_time / num_reports;

                printf("Average node to base communication time: %f seconds \n", avg_base_comm_time);

                FILE *file = fopen("log.txt", "a"); // "a" means append, creates the file if it doesn't exist

                fprintf(file, "Total node to node communication time: %f seconds \n", node_total_comm_time);
                fprintf(file, "Average node to node communication time: %f seconds \n", avg_node_comm_time);
                fprintf(file, "Total node to base communication time: %f seconds \n", total_comm_time);
                fprintf(file, "Average node to base communication time: %f seconds \n", avg_base_comm_time);

                fclose(file);
                
                terminate = 1;

            }
        }
        
    }

    return 0;
}

int runChargingNode(MPI_Comm world_comm, MPI_Comm split_comm, int dims[2], int NUM_NODES) {
    int my_rank, size; 
    int terminate = 0; // boolean to end execution
    int ierr; // to check if grid produced correctly
    int periods[2] = {0, 0}; // periods array to make sure it doesnt wrap
    MPI_Comm grid_comm; // new communicator for grid
    int coord[2];
    int nbr_down, nbr_up, nbr_left, nbr_right;
    int ports[NUM_PORTS] = {0};
    Table table = {.count = 0}; // Initialize 

    int report_flag = 0; // checks if report needs to be sent


    ierr = MPI_Cart_create(split_comm, 2, dims, periods, 0, &grid_comm);
    if (ierr != 0) printf("ERROR[%d] creating CART\n", ierr);

    MPI_Comm_rank(grid_comm, &my_rank);
    MPI_Comm_size(grid_comm, &size);

    MPI_Cart_coords(grid_comm, my_rank, 2, coord);
    MPI_Cart_shift(grid_comm, 0, 1, &nbr_up, &nbr_down); // gets top and bottom neighbours
    MPI_Cart_shift(grid_comm, 1, 1, &nbr_left, &nbr_right); // gets left and right neighbours

    int neighbours[4] = {nbr_up, nbr_down, nbr_left, nbr_right};

    int my_availability = checkNodeAvailability(ports); 

    int nbr_availabilities[4] = {-1, -1, -1, -1};

    MPI_Datatype mpi_report_type = create_mpi_report();

    double total_comm_time = 0;
    double comm_time, start_time, end_time;
    int num_communications = 0;


    #pragma omp parallel num_threads(NUM_PORTS + 3), shared(ports, table, my_availability, report_flag, nbr_availabilities, terminate, total_comm_time, num_communications)
    {
        int thread_id = omp_get_thread_num();
        while (terminate == 0) {
            if (thread_id >= 0 && thread_id <= NUM_PORTS - 1) {
                #pragma omp critical 
                {    
                    ports[thread_id] = rand() % 2;
                }
            }

            my_availability = checkNodeAvailability(ports); 

            if (thread_id == 0) {

                time_t t;
                struct tm tm;
                time(&t);
                tm = *localtime(&t);

                MPI_Request send_request[4]; // array for send requests
                MPI_Request receive_request[4]; // array for receive requests
                MPI_Status send_status[4]; // array for send statuses
                MPI_Status receive_status[4]; // array for receive statuses

                int num_requests = 0; // tracks number of requests sent out

                TableEntry entry;
                
                entry.year = tm.tm_year + 1900;
                entry.month = tm.tm_mon + 1;
                entry.day = tm.tm_mday;
                entry.hour = tm.tm_hour;
                entry.minute = tm.tm_min;
                entry.second = tm.tm_sec;
                entry.availability = my_availability;

                insertTable(&table, entry);
                
                int full_flag = 0; 

                // 1d) checks if full
                if (my_availability <= FULL_THRESHOLD) {
                    
                    full_flag = 1;

                    start_time = MPI_Wtime();

                    for (int i = 0; i < 4; i++) {
                        if (neighbours[i] != -2) {
                            MPI_Isend(&my_rank, 1, MPI_INT, neighbours[i], REQUEST_TAG, grid_comm, &send_request[num_requests]); // send request
                            MPI_Irecv(&nbr_availabilities[i], 1, MPI_INT, neighbours[i], AVAILABILITY_TAG, grid_comm, &receive_request[num_requests]); // receive availability
                            num_requests++;
                        }
                    }

                    MPI_Waitall(num_requests, send_request, send_status);
                    MPI_Waitall(num_requests, receive_request, receive_status);

                    end_time = MPI_Wtime();    

                    comm_time = end_time - start_time;

                    total_comm_time += comm_time;

                    num_communications++;

                }

                if (full_flag == 1) {

                    report_flag = 0; // signals to send report when == 1
                    int all_full = 1; // variable to check if all neighbours are full

                    for (int i = 0; i < 4; i++) {  // Assuming there are 4 neighbors
                        if (nbr_availabilities[i] > FULL_THRESHOLD) {
                            all_full = 0; // all neighbours are not full
                            // commented out for testing
                            // printf("Neighbour %d has availability %d\n", neighbours[i], nbr_availabilities[i]);
                        }
                    }
                    if (all_full) {
                        printf("neighbour availabilities %d : %d, %d, %d, %d\n", my_rank, nbr_availabilities[0], nbr_availabilities[1], nbr_availabilities[2], nbr_availabilities[3]);
                        #pragma omp critical
                        {
                            report_flag = 1;
                        }
                    }
                    
                    
                } 

                sleep(CYCLE_TIME);

                

            } else if (thread_id == NUM_PORTS) {
                while (terminate == 0) {
                    int flag = 0;
                    int req_rank;
                    MPI_Iprobe(MPI_ANY_SOURCE, REQUEST_TAG, grid_comm, &flag, MPI_STATUS_IGNORE);
                    if (flag) {
                        MPI_Recv(&req_rank, 1, MPI_INT, MPI_ANY_SOURCE, REQUEST_TAG, grid_comm, MPI_STATUS_IGNORE);
                        MPI_Send(&my_availability, 1, MPI_INT, req_rank, AVAILABILITY_TAG, grid_comm);
                    }
                }

            } else if (thread_id == NUM_PORTS + 1) {
                int num_neighbours = 0;
                for (int i = 0; i < 4; i++) {
                    if (neighbours[i] != -2) {
                        num_neighbours++;
                    }
                }

                while (terminate == 0) {
                    if (report_flag == 1) {
                        // Reset report flag
                        report_flag = 0;

                        // Create and populate the report
                        Report report;

                        // Populate the current date and time
                        time_t t;
                        time(&t);
                        report.year = table.data[0].year;
                        report.month = table.data[0].month;
                        report.day = table.data[0].day;
                        report.hour = table.data[0].hour;
                        report.minute = table.data[0].minute;
                        report.second = table.data[0].second;

                        // Add start time
                        report.start_time = MPI_Wtime();

                        // Populate the coordinates
                        report.coords[0] = coord[0];
                        report.coords[1] = coord[1];

                        // populate availability
                        report.availability = table.data[0].availability;

                        // Populate the neighbours
                        for (int i = 0; i < 4; i++) {
                            report.neighbours[i] = neighbours[i];
                        }

                        // Populate the neighbour availabilities
                        for (int i = 0; i < 4; i++) {
                            report.nbr_availabilities[i] = nbr_availabilities[i];
                        }

                        // populate the number of messages exchanged
                        report.numMessagesExchanged = num_neighbours * 2;

                        #pragma omp critical
                        {
                            // Send the report (assuming destination rank is known and tag is defined)
                            MPI_Send(&report, 1, mpi_report_type, 0, REPORT_TAG, MPI_COMM_WORLD); //base station is 0 in MPI_COMM_WORLD
                            int* nearby_nodes = (int*)malloc(NUM_NODES * sizeof(int));
                            MPI_Recv(nearby_nodes, NUM_NODES, MPI_INT, 0, NEARBY_NODES_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                            int no_nearby = 1;

                            for (int i = 0; i < NUM_NODES; i++) {
                                if (nearby_nodes[i] == 1) {
                                    no_nearby = 0;
                                    printf("cars attempting to charge at node %d can instead charge at node %d\n",my_rank, i);
                                }
                            }
                            if (no_nearby == 1) {
                                    printf("There are no available nodes near %d.\n", my_rank);
                                }
                        }
                    }
                }
            } else if (thread_id == NUM_PORTS + 2) {
                int terminate_flag = 0;
                while (terminate == 0) {
                    MPI_Iprobe(0, TERMINATION_TAG, MPI_COMM_WORLD, &terminate_flag, MPI_STATUS_IGNORE);
                    if (terminate_flag) {
                        
                        CommTime comm_time;
                        comm_time.total_comm_time = total_comm_time;
                        comm_time.num_communications = num_communications;

                        printf("num communications: %d\n", comm_time.num_communications);

                        MPI_Datatype mpi_commtime_type = create_mpi_commtime();

                        MPI_Send(&comm_time, 1, mpi_commtime_type, 0, COMM_TIME_TAG, MPI_COMM_WORLD);

                        MPI_Type_free(&mpi_commtime_type);  


                        terminate = 1;
                    }
                }
            }
                
        }
    }
    printf("%d terminate\n", my_rank);
    MPI_Type_free(&mpi_report_type);
    return 0;
}

int main(int argc, char* argv[]) {

    int world_rank, world_size, nrow, ncol;
    int dims[2]; // dimensions array

    MPI_Comm split_comm; // communicator for split of base station and nodes
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);


    if (argc == 3) {
		nrow = atoi (argv[1]);
		ncol = atoi (argv[2]);
		dims[0] = nrow; /* number of rows */
		dims[1] = ncol; /* number of columns */
	} else {
		nrow=ncol=(int)sqrt(world_size);
		dims[0]=dims[1]=0;
	}

    FILE *file = fopen("log.txt", "w");
    if (file == NULL) {
        perror("Error opening file");
        return -1;
    }
    fclose(file);

    int isBaseStation = (world_rank == 0);
    MPI_Comm_split(MPI_COMM_WORLD, isBaseStation, 0, &split_comm); // split base station and charging nodes

    int NUM_NODES = nrow*ncol;

    if (isBaseStation) {
        runBaseStation(NUM_NODES);
    } else {
        runChargingNode(MPI_COMM_WORLD, split_comm, dims, NUM_NODES);
    }

    MPI_Finalize();

    printf("FINISH\n");
    

    return 0;
}
