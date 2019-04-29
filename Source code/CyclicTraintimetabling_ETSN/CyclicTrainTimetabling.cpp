//cyclic train timetabling based on extended time-discritized space-time network

#include "stdafx.h"
#include <iostream>
#include <fstream>
#include <list> 
#include <omp.h>
#include <algorithm>
#include <time.h>
#include "CSVParser.h"
#include <functional>
#include <stdio.h>   
#include <tchar.h>
#include <windows.h>
#include <vector>
#include <random>
#include <chrono>

#define _MAX_LABEL_COST 999999
#define _MAX_STATES 1

#define _MAX_NUMBER_OF_TIME_INTERVALS 720

// Linear congruential generator 
#define LCG_a 17364
#define LCG_c 0
#define LCG_M 65521


using namespace std;

TCHAR g_SettingFileName[_MAX_PATH] = _T("./Settings.txt");

FILE* g_pFileDebugLog = NULL;
FILE* g_pFileDebugLog_LR = NULL;
FILE* g_pFileDebugLog_ADMM = NULL;
FILE* g_pFileOutputLog = NULL;
FILE* g_LR_iteration_Log = NULL;
FILE* g_ADMM_iteration_Log = NULL;
FILE* g_LR_algorithmic_times_Log = NULL;
FILE* g_ADMM_algorithmic_times_Log = NULL;

int g_cycle_length = 120;
int g_number_of_simulation_intervals = 720;
int g_number_of_intervals_in_master_schedule = 360;

CTime g_SolutionStartTime;
clock_t LR_Initialization_start, LR_Initialization_end;
clock_t LR_LB_start, LR_LB_end;
clock_t LR_UB_start, LR_UB_end;
clock_t LR_LMU_start, LR_LMU_end;

float LR_Initialization_time = 0;
float LR_LB_time = 0;
float LR_UB_time = 0;
float LR_LMU_time = 0;

clock_t ADMM_Initialization_start, ADMM_Initialization_end;
clock_t ADMM_LB_start, ADMM_LB_end;
clock_t ADMM_UB_start, ADMM_UB_end;
clock_t ADMM_LMU_start, ADMM_LMU_end;
float ADMM_Initialization_time = 0;
float ADMM_LB_time = 0;
float ADMM_UB_time = 0;
float ADMM_LMU_time = 0;

int g_dp_algorithm_debug_flag = 0;
int g_LR_algorithm_debug_flag = 3; //default 3
int g_ADMM_algorithm_debug_flag = 3;
int g_upper_bound_solution_check_flag = 1;
int g_deduce_feasible_upperbound_in_ADMM_flag = 0;
bool g_output_log_flag = true;
int H = 3; // the number of copies

//headway requirements
int g_departure_headway_stop = 5;
int g_departure_headway_passing = 3;
int g_arrival_headway_stop = 4;
int g_arrival_headway_passing = 3;

int g_accelerating_time_fast_train = 3; //speed grade: 1
int g_decelerating_time_fast_train = 3; //speed grade: 1
int g_accelerating_time_slow_train = 2; //speed grade: 2
int g_decelerating_time_slow_train = 3; //speed grade: 2

int g_number_of_links = 0;
int g_number_of_nodes = 0;
int g_number_of_zones = 0;

int g_number_of_LR_iterations = 0;
int g_number_of_ADMM_iterations = 1000;

int penalty_parameter_increasing_strategy = 1; //0, simple; 1, complex
int number_of_lagrangian_iterations_for_ADMM = 200; //used in strategy 1
int g_penalty_gamma_ratio = 2;
vector<int> g_penalty_term_vector;

float g_penalty_RHO_initial = 4;
float g_penalty_RHO_incremental = 2;
float g_gama_ratio = 1;
int g_penalty_RHO_update_iteration = 20;
int g_ADMM_feasible_and_stop_flag = 1;
bool g_reoptimization_flag = false;
bool g_random_permutation = false;

float g_best_upper_bound = 99999;
float g_best_lower_bound = -99999;
int g_primal_residual_error = 99999;

float g_best_ADMM_Feasible_lower_bound = 99999;
float optimality_gap = 0;
int freStop_best_feas_Upper_bound_flag = 1;
int freStop_best_feas_Lower_bound_flag = 0;
float g_stepSize = 1;
float g_minimum_subgradient_step_size = 0.1;
float g_initial_LR_multiplier_value = 2;

int siding_track_type = 11;
int dummy_track_type = 12;

int g_number_of_agents;
int g_number_of_threads = 1;

float g_number_of_seconds_per_interval = 0.2;
int g_number_of_intervals_per_min = 60 / g_number_of_seconds_per_interval;
int g_number_of_simulation_minutes = 100;

vector <vector<int>> g_candidate_frequence_stop_pattern;
vector<float> g_LR_FreStop_best_lower_bound;
vector<float> g_LR_FreStop_best_Upper_bound;
vector<float> g_LR_FreStop_best_Optimality_gap;
vector<int> g_LR_FreStop_best_feas_Upper_bound_flag;

vector<float> g_ADMM_FreStop_best_lower_bound;
vector<float> g_ADMM_FreStop_best_Upper_bound;
vector<float> g_ADMM_FreStop_best_Optimality_gap;
vector<int> g_ADMM_FreStop_best_feas_Lower_bound_flag;

vector<int> g_agent_sequence_by_LP_ADMM; //train sequences based on lagrangian profit
vector<int> g_ADMM_optimization_sequence; //the global optimization sequences in ADMM

int m_ListFront;
int m_ListTail;

std::map<int, int> g_link_key_to_seq_no_map;  // hush table, map key to internal link sequence no.
std::map<int, int> g_internal_link_no_map;  // hash table
std::map<int, int> g_external_link_id_map;  // hash table

std::map<int, int> g_internal_agent_no_map; // map exteranl agent id to internal agent no.
std::map<int, int> g_external_agent_id_map; // map internal agent no. to external agent id

std::map<int, int> g_internal_node_seq_no_map;  // hush table, map external node number to internal node sequence no. 
std::map<int, int> g_internal_node_seq_no_to_node_id_map;  // hush table, map external node number to internal node sequence no. 

std::map<int, int> g_internal_signal_seq_no_map;
std::map<int, int> g_internal_signal_seq_no_to_node_id_map;

template <typename T>
T **Allocate2DDynamicArray(int nRows, int nCols)
{
	T **dynamicArray;

	dynamicArray = new T*[nRows];

	for (int i = 0; i < nRows; i++)
	{
		dynamicArray[i] = new T[nCols];

		if (dynamicArray[i] == NULL)
		{
			cout << "Error: insufficent memory.";
			g_ProgramStop();
		}
	}

	return dynamicArray;
}

template <typename T>
void Deallocate2DDynamicArray(T** dArray, int nRows)
{
	for (int x = 0; x < nRows; x++)
	{
		delete[] dArray[x];
	}

	delete[] dArray;
}

template <typename T>
T ***Allocate3DDynamicArray(int nX, int nY, int nZ)
{
	T ***dynamicArray;

	dynamicArray = new (std::nothrow) T**[nX];

	if (dynamicArray == NULL)
	{
		cout << "Error: insufficient memory.";
		g_ProgramStop();
	}

	for (int x = 0; x < nX; x++)
	{
		dynamicArray[x] = new (std::nothrow) T*[nY];

		if (dynamicArray[x] == NULL)
		{
			cout << "Error: insufficient memory.";
			g_ProgramStop();
		}

		for (int y = 0; y < nY; y++)
		{
			dynamicArray[x][y] = new (std::nothrow) T[nZ];
			if (dynamicArray[x][y] == NULL)
			{
				cout << "Error: insufficient memory.";
				g_ProgramStop();
			}
		}
	}

	return dynamicArray;
}

template <typename T>
void Deallocate3DDynamicArray(T*** dArray, int nX, int nY)
{
	if (!dArray)
		return;
	for (int x = 0; x < nX; x++)
	{
		for (int y = 0; y < nY; y++)
		{
			delete[] dArray[x][y];
		}

		delete[] dArray[x];
	}

	delete[] dArray;
}

long g_GetLinkSeqNo(int from_node_no, int to_node_no)
{
	long link_key = from_node_no * 100000 + to_node_no;

	if (g_link_key_to_seq_no_map.find(link_key) != g_link_key_to_seq_no_map.end())
		return g_link_key_to_seq_no_map[link_key];
	else
		return -1;
}

class CLink
{
public:
	CLink()
	{
		cost = 0;
		free_flow_travel_time_in_min = 1;
		speed_limit = 10;
	}

	int** departure_state_dependent_travel_time_matrix;
	int** arrival_state_dependent_travel_time_matrix;

	float** state_time_dependent_departure_LR_multiplier_matrix;
	float** state_time_dependent_arrival_LR_multiplier_matrix;

	int** state_dependent_time_dependent_lastiternum;

	//ADMM_multiplier_matrix is used in the searching process and updated by LR_multiplier_matrix
	float** state_dependent_time_dependent_departure_ADMM_multiplier_matrix;
	float** state_dependent_time_dependent_arrival_ADMM_multiplier_matrix;

	int* time_depedent_capacity_matrix;

	int* g_link_time_departure_visit_counts; //(i, j, t)
	int* g_link_time_arrival_visit_counts; //(i, j, t)

	int** g_link_time_departure_train_visit_flag; //(i, j, t, a)
	int** g_link_time_arrival_train_visit_flag; //(i, j, t, a)

	int* g_link_time_visit_counts_departure_in_upper_bound; //(i, j, t)
	int* g_link_time_visit_counts_arrival_in_upper_bound; //(i, j, t)

	int** g_link_time_train_visit_flag_departure_in_upper_bound; //(i, j, t, a)
	int** g_link_time_train_visit_flag_arrival_in_upper_bound; //(i, j, t, a)

	void Setup_State_Dependent_Data_Matrix()
	{
		departure_state_dependent_travel_time_matrix = Allocate2DDynamicArray<int>(g_number_of_simulation_intervals, _MAX_STATES);
		arrival_state_dependent_travel_time_matrix = Allocate2DDynamicArray<int>(g_number_of_simulation_intervals, _MAX_STATES);

		time_depedent_capacity_matrix = new int[g_number_of_simulation_intervals];

		g_link_time_departure_train_visit_flag = Allocate2DDynamicArray<int>(g_number_of_simulation_intervals, g_number_of_agents);
		g_link_time_arrival_train_visit_flag = Allocate2DDynamicArray<int>(g_number_of_simulation_intervals, g_number_of_agents);

		g_link_time_departure_visit_counts = new int[g_number_of_simulation_intervals];
		g_link_time_arrival_visit_counts = new int[g_number_of_simulation_intervals];

		g_link_time_train_visit_flag_departure_in_upper_bound = Allocate2DDynamicArray<int>(g_number_of_simulation_intervals, g_number_of_agents);
		g_link_time_train_visit_flag_arrival_in_upper_bound = Allocate2DDynamicArray<int>(g_number_of_simulation_intervals, g_number_of_agents);

		g_link_time_visit_counts_departure_in_upper_bound = new int[g_number_of_simulation_intervals];
		g_link_time_visit_counts_arrival_in_upper_bound = new int[g_number_of_simulation_intervals];

		for (int t = 0; t < g_number_of_simulation_intervals; t++)
		{
			if (link_type == siding_track_type || link_type == dummy_track_type)
			{
				time_depedent_capacity_matrix[t] = station_track_capacity;
			}
			else
			{
				time_depedent_capacity_matrix[t] = 1;
			}

			g_link_time_departure_visit_counts[t] = 0;
			g_link_time_arrival_visit_counts[t] = 0;

			g_link_time_visit_counts_departure_in_upper_bound[t] = 0;
			g_link_time_visit_counts_arrival_in_upper_bound[t] = 0;

			for (int s = 0; s < _MAX_STATES; s++)
			{
				if (link_type == siding_track_type || link_type == dummy_track_type)
				{
					departure_state_dependent_travel_time_matrix[t][s] = station_track_capacity - 1;
					arrival_state_dependent_travel_time_matrix[t][s] = station_track_capacity - 1;
				}
				else
				{
					departure_state_dependent_travel_time_matrix[t][s] = 0;
					arrival_state_dependent_travel_time_matrix[t][s] = 0;
				}
			}

			for (int a = 0; a < g_number_of_agents; a++)
			{
				g_link_time_departure_train_visit_flag[t][a] = 0;
				g_link_time_arrival_train_visit_flag[t][a] = 0;

				g_link_time_train_visit_flag_departure_in_upper_bound[t][a] = 0;
				g_link_time_train_visit_flag_arrival_in_upper_bound[t][a] = 0;
			}
		}

		state_time_dependent_departure_LR_multiplier_matrix = Allocate2DDynamicArray<float>(g_number_of_simulation_intervals, _MAX_STATES);
		state_time_dependent_arrival_LR_multiplier_matrix = Allocate2DDynamicArray<float>(g_number_of_simulation_intervals, _MAX_STATES);

		state_dependent_time_dependent_lastiternum = Allocate2DDynamicArray<int>(g_number_of_simulation_intervals, _MAX_STATES);

		state_dependent_time_dependent_departure_ADMM_multiplier_matrix = Allocate2DDynamicArray<float>(g_number_of_simulation_intervals, _MAX_STATES);
		state_dependent_time_dependent_arrival_ADMM_multiplier_matrix = Allocate2DDynamicArray<float>(g_number_of_simulation_intervals, _MAX_STATES);

		for (int t = 0; t < g_number_of_simulation_intervals; t++)
		{
			for (int s = 0; s < _MAX_STATES; s++)
			{
				state_time_dependent_departure_LR_multiplier_matrix[t][s] = g_initial_LR_multiplier_value;
				state_time_dependent_arrival_LR_multiplier_matrix[t][s] = g_initial_LR_multiplier_value;

				state_dependent_time_dependent_lastiternum[t][s] = 0;

				state_dependent_time_dependent_departure_ADMM_multiplier_matrix[t][s] = 0;
				state_dependent_time_dependent_arrival_ADMM_multiplier_matrix[t][s] = 0;
			}
		}
	}

	int link_id;
	int same_link_id;
	int station_track_capacity;

	int link_type;
	int lane_id;
	string name;

	int shp_id;
	int from_node_id; // external node numbers, 
	int to_node_id;
	int direction; // 1 is upward direction; 2 is downward direction

	int link_seq_no;
	int from_node_seq_no; // starting from 0, sequential numbers 
	int to_node_seq_no;
	float cost;

	float free_flow_travel_time_in_min;
	float length;
	float speed_limit;
	int type;
	float travel_time;

	double x;
	double y;
	double local_y;
};

class CNode
{
public:

	void Setup_Node_State_Matrix()
	{
		node_state_matrix = Allocate2DDynamicArray<int>(g_number_of_simulation_intervals, _MAX_STATES);

		for (int t = 0; t < g_number_of_simulation_intervals; t++)
		{
			for (int s = 0; s < _MAX_STATES; s++)
			{
				node_state_matrix[t][s] = 1;
			}
		}
	}

	int node_seq_no;  // sequence number
	int node_id;  //external node number
	string name;
	double x;
	double y;

	int station_id;
	int station_sequence;

	int origin_destination_flag;

	int** node_state_matrix;

	std::vector<CLink> m_outgoing_node_vector;
};

std::vector<CNode> g_node_vector;
std::vector<CLink> g_link_vector;

class CAgent
{
public:
	unsigned int m_RandomSeed;

	CAgent()
	{
		agent_vector_seq_no = -1;
		path_index = -1;
	}

	float GetRandomRatio()
	{
		m_RandomSeed = (LCG_a * m_RandomSeed + LCG_c) % LCG_M;  //m_RandomSeed is automatically updated.

		return float(m_RandomSeed) / LCG_M;
	}

	int agent_id;
	int agent_vector_seq_no;
	int origin_node_id;
	int destination_node_id;
	int path_index;

	float earliest_departure_time;
	int departure_time_in_simulation_interval;
	float departure_time_window;
	int frequency;
	int min_frequency;
	int max_frequency;

	float free_flow_travel_time;
	float travel_time_in_dual_solution;
	float speed_grade;

	int direction;

	vector<int> set_of_allowed_links_LR;
	vector<int> m_set_of_allowed_links_flag_LR;

	float travel_time_in_min;
	vector<int> partial_schedule_node_vector;
	vector<int> partial_schedule_node_TA; // -1 means non-existence
	vector<int> partial_schedule_node_TD; // -1 means non-existence

	std::vector<int> path_link_seq_no_vector;
	std::vector<int> path_node_id_vector;
	std::vector<int> path_timestamp_vector;

	std::vector<int> path_new_node_id_vector;
	std::vector<int> path_new_link_id_vector;
	std::vector<int> path_new_timestamp_vector;
	std::vector<int> path_link_TA_vector;
	std::vector<int> path_link_TD_vector;

	std::vector<int> path_new_link_id_vector_temp;
	std::vector<int> path_link_TA_vector_temp;
	std::vector<int> path_link_TD_vector_temp;

	std::vector<int> g_output_link_NO_vector;
	std::vector<int> g_output_link_TA_vector;
	std::vector<int> g_output_link_TD_vector;

	std::vector<int> path_link_seq_no_vector_upper_bound;
	std::vector<int> path_timestamp_vector_upper_bound;
	std::vector<int> path_node_id_vector_upper_bound;
	std::vector<int> path_new_link_id_vector_upper_bound;
	std::vector<int> path_link_TA_vector_upper_bound;
	std::vector<int> path_link_TD_vector_upper_bound;
	std::vector<int> path_new_node_id_vector_upper_bound;
	std::vector<int> path_new_timestamp_vector_upper_bound;

	int* min_link_travel_time;
	int* max_link_travel_time;
	int m_current_link_seq_no;
	float m_DeviationRatio;

	void Setup_Time_Dependent_Data_Matrix()
	{
		min_link_travel_time = new int[g_number_of_links];
		max_link_travel_time = new int[g_number_of_links];

		for (int l = 0; l < g_number_of_links; l++)
		{
			min_link_travel_time[l] = 1;
			max_link_travel_time[l] = 1;
		}
	}

	bool operator<(const CAgent &other) const
	{
		return earliest_departure_time < other.earliest_departure_time;
	}

};

vector<CAgent> g_agent_vector;
std::map<int, int> g_map_agent_id_to_agent_vector_seq_no;

void g_ProgramStop()
{
	cout << "TimetableLite Program stops. Press any key to terminate. Thanks!" << endl;
	getchar();
	exit(0);
};

//split the string by ";"
vector<string> split(const string &s, const string &seperator) {
	vector<string> result;
	typedef string::size_type string_size;
	string_size i = 0;

	while (i != s.size()) {
		int flag = 0;
		while (i != s.size() && flag == 0) {
			flag = 1;
			for (string_size x = 0; x < seperator.size(); ++x)
				if (s[i] == seperator[x]) {
					++i;
					flag = 0;
					break;
				}
		}

		flag = 0;
		string_size j = i;
		while (j != s.size() && flag == 0) {
			for (string_size x = 0; x < seperator.size(); ++x)
				if (s[j] == seperator[x]) {
					flag = 1;
					break;
				}
			if (flag == 0)
				++j;
		}
		if (i != j) {
			result.push_back(s.substr(i, j - i));
			i = j;
		}
	}
	return result;
}

void g_ReadInputData()
{
	// initialize  the counter to 0
	g_number_of_nodes = 0;
	g_number_of_links = 0;

	// step 1: read node file 
	CCSVParser parser;

	if (parser.OpenCSVFile("input_node.csv", true))
	{
		int internal_node_seq_no = 0;
		double x, y;

		std::map<int, int> node_id_map;

		while (parser.ReadRecord()) // if this line contains [] mark, then we will also read field headers.
		{
			string name;

			int node_type;
			int node_id;

			if (parser.GetValueByFieldName("node_id", node_id) == false)
				continue;

			if (g_internal_node_seq_no_map.find(node_id) != g_internal_node_seq_no_map.end())
			{
				continue; //has been defined
			}

			g_internal_node_seq_no_map[node_id] = internal_node_seq_no;
			g_internal_node_seq_no_to_node_id_map[internal_node_seq_no] = node_id;

			parser.GetValueByFieldName("x", x, false);
			parser.GetValueByFieldName("y", y, false);
			parser.GetValueByFieldName("name", name);

			CNode node;  // create a node object
			node.node_id = node_id;
			node.node_seq_no = internal_node_seq_no;
			node.name = name;

			parser.GetValueByFieldName("station_id", node.station_id, false);
			parser.GetValueByFieldName("station_sequence", node.station_sequence, false);

			node.x = x;
			node.y = y;
			internal_node_seq_no++;

			g_node_vector.push_back(node);  // push it to the global node vector

			g_number_of_nodes++;
			if (g_number_of_nodes % 1000 == 0)
				cout << "reading " << g_number_of_nodes << " nodes.. " << endl;
		}

		cout << "number of nodes = " << g_number_of_nodes << endl;

		fprintf(g_pFileOutputLog, "number of nodes =,%d\n", g_number_of_nodes);
		parser.CloseCSVFile();
	}

	// step 2: read link file 
	CCSVParser parser_link;

	if (parser_link.OpenCSVFile("input_link.csv", true))
	{
		while (parser_link.ReadRecord())  // if this line contains [] mark, then we will also read field headers.
		{
			CLink link;  // create a link object

			if (parser_link.GetValueByFieldName("link_id", link.link_id) == false)
				continue;
			if (parser_link.GetValueByFieldName("from_node_id", link.from_node_id) == false)
				continue;
			if (parser_link.GetValueByFieldName("to_node_id", link.to_node_id) == false)
				continue;

			parser_link.GetValueByFieldName("name", link.name);
			parser_link.GetValueByFieldName("direction", link.direction);
			parser_link.GetValueByFieldName("link_type", link.link_type);
			parser_link.GetValueByFieldName("lane_id", link.lane_id);
			parser_link.GetValueByFieldName("same_link_id", link.same_link_id);
			parser_link.GetValueByFieldName("station_track_capacity", link.station_track_capacity);

			// add the to node id into the outbound (adjacent) node list
			int internal_from_node_seq_no = g_internal_node_seq_no_map[link.from_node_id];  // map external node number to internal node seq no. 
			int internal_to_node_seq_no = g_internal_node_seq_no_map[link.to_node_id];

			link.from_node_seq_no = internal_from_node_seq_no;
			link.to_node_seq_no = internal_to_node_seq_no;
			link.link_seq_no = g_number_of_links;

			parser_link.GetValueByFieldName("local_y", link.local_y);

			g_internal_link_no_map[link.link_id] = link.link_seq_no;
			g_external_link_id_map[link.link_seq_no] = link.link_id;

			float length = 1; // km or mile

			parser_link.GetValueByFieldName("length", length);
			parser_link.GetValueByFieldName("speed_limit", link.speed_limit);

			link.free_flow_travel_time_in_min = 60 * length / link.speed_limit;

			link.length = length;
			link.cost = length / link.speed_limit * 60; //calculate link cost based length and speed limit

			link.x = (g_node_vector[link.from_node_seq_no].x + g_node_vector[link.to_node_seq_no].x) / 2.0;
			link.y = (g_node_vector[link.from_node_seq_no].y + g_node_vector[link.to_node_seq_no].y) / 2.0;

			g_node_vector[internal_from_node_seq_no].m_outgoing_node_vector.push_back(link);  // add this link to the corresponding node as part of outgoing node/link

			long link_key = internal_from_node_seq_no * 100000 + internal_to_node_seq_no;

			g_link_key_to_seq_no_map[link_key] = link.link_seq_no;

			g_link_vector.push_back(link);

			g_number_of_links++;

			if (g_number_of_links % 1000 == 0)
				cout << "reading " << g_number_of_links << " links.. " << endl;
		}
	}

	cout << "number of links = " << g_number_of_links << endl;

	fprintf(g_pFileOutputLog, "number of links =,%d\n", g_number_of_links);

	parser_link.CloseCSVFile();

	g_number_of_agents = 0;

	if (g_number_of_agents == 0)
	{
		CCSVParser parser_agent;
		std::vector<int> path_node_sequence;
		string path_node_sequence_str;

		std::vector<int> path_schedule_time_sequence;
		string path_schedule_time_sequence_str;

		if (parser_agent.OpenCSVFile("input_agent.csv", true))
		{
			while (parser_agent.ReadRecord())
			{
				CAgent agent;  // create an agent object
				if (parser_agent.GetValueByFieldName("agent_id", agent.agent_id) == false)
					continue;

				agent.m_RandomSeed = agent.agent_id;

				g_internal_agent_no_map[agent.agent_id] = g_number_of_agents;
				g_external_agent_id_map[g_number_of_agents] = agent.agent_id;
				agent.agent_vector_seq_no = g_internal_agent_no_map[agent.agent_id];

				int origin_node_id = 0;
				int destination_node_id = 0;

				parser_agent.GetValueByFieldName("from_origin_node_id", origin_node_id, false);
				agent.origin_node_id = origin_node_id;

				parser_agent.GetValueByFieldName("to_destination_node_id", destination_node_id, false);
				agent.destination_node_id = destination_node_id;

				if (g_internal_node_seq_no_map.find(origin_node_id) == g_internal_node_seq_no_map.end() || g_internal_node_seq_no_map.find(destination_node_id) == g_internal_node_seq_no_map.end())
					continue;

				parser_agent.GetValueByFieldName("earliest_departure_time", agent.earliest_departure_time);

				parser_agent.GetValueByFieldName("frequency", agent.frequency);
				agent.departure_time_window = ceil(g_cycle_length / agent.frequency);

				parser_agent.GetValueByFieldName("direction", agent.direction);

				//speed grade: 1 is fast train; 2 is slow train
				parser_agent.GetValueByFieldName("speed_grade", agent.speed_grade);

				agent.departure_time_in_simulation_interval = agent.earliest_departure_time;

				string set_of_allowed_links_str;
				vector<string> set_of_allowed_links_sub_str;
				parser_agent.GetValueByFieldName("set_of_allowed_links", set_of_allowed_links_str);

				//initialize the travel time and cost matrix
				agent.Setup_Time_Dependent_Data_Matrix();

				for (int l = 0; l < g_number_of_links; l++)
				{
					agent.m_set_of_allowed_links_flag_LR.push_back(0);
				}

				int temp_link_no = 0;

				//set_of_allowed_links
				if (set_of_allowed_links_str != "")
				{
					set_of_allowed_links_sub_str = split(set_of_allowed_links_str, ";");
					int set_of_allowed_links_size = max(set_of_allowed_links_sub_str.size(), 1);

					for (int i = 0; i < set_of_allowed_links_size; i++)
					{
						temp_link_no = g_internal_link_no_map[stoi(set_of_allowed_links_sub_str[i])]; //stoi: string to int

						agent.set_of_allowed_links_LR.push_back(temp_link_no);
						agent.m_set_of_allowed_links_flag_LR[temp_link_no] = 1;
					}
				}

				string min_link_travel_time_str;
				vector<string> min_link_travel_time_sub_str;
				parser_agent.GetValueByFieldName("min_link_travel_time", min_link_travel_time_str);

				int temp_time = 1;

				//min_link_travel_time
				if (min_link_travel_time_str != "")
				{
					min_link_travel_time_sub_str = split(min_link_travel_time_str, ";");
					int min_link_travel_time_size = max(min_link_travel_time_sub_str.size(), 1);

					for (int i = 0; i < g_link_vector.size(); i++)
					{
						temp_time = stoi(min_link_travel_time_sub_str[i]); //stoi: string to int
						temp_link_no = g_link_vector[i].link_seq_no;
						agent.min_link_travel_time[temp_link_no] = temp_time;
					}
				}

				//max_link_travel_time
				string max_link_travel_time_str;
				vector<string> max_link_travel_time_sub_str;
				parser_agent.GetValueByFieldName("max_link_travel_time", max_link_travel_time_str);

				if (max_link_travel_time_str != "")
				{
					max_link_travel_time_sub_str = split(max_link_travel_time_str, ";");
					int max_link_travel_time_size = max(max_link_travel_time_sub_str.size(), 1);

					for (int i = 0; i < g_link_vector.size(); i++)
					{
						temp_time = stoi(max_link_travel_time_sub_str[i]); //stoi: string to int
						temp_link_no = g_link_vector[i].link_seq_no;
						agent.max_link_travel_time[temp_link_no] = temp_time;
					}
				}

				int next_link_no = 0;
				int former_link_no = 0;

				//revise the min and max travel time by acceleration/deceleration time
				for (int l = 0; l < agent.set_of_allowed_links_LR.size(); l++)  // for each link in the path of this agent
				{
					temp_link_no = agent.set_of_allowed_links_LR[l];

					if (g_link_vector[temp_link_no].link_type != siding_track_type && g_link_vector[temp_link_no].link_type != dummy_track_type) // 11 is the station track
					{
						if (l == 0) //first section
						{
							//acceleration time
							if (agent.speed_grade == 1) //'G', fast train
							{
								agent.min_link_travel_time[temp_link_no] += g_accelerating_time_fast_train;
								agent.max_link_travel_time[temp_link_no] += g_accelerating_time_fast_train;
							}
							else if (agent.speed_grade == 2) //'D', slow train
							{
								agent.min_link_travel_time[temp_link_no] += g_accelerating_time_slow_train;
								agent.max_link_travel_time[temp_link_no] += g_accelerating_time_slow_train;
							}

							//deceleration time
							next_link_no = agent.set_of_allowed_links_LR[l + 1];

							if (g_link_vector[next_link_no].link_type == siding_track_type) //train will stop
							{
								if (agent.speed_grade == 1) //'G', fast train
								{
									agent.min_link_travel_time[temp_link_no] += g_decelerating_time_fast_train;
									agent.max_link_travel_time[temp_link_no] += g_decelerating_time_fast_train;
								}
								else if (agent.speed_grade == 2) //'D', slow train
								{
									agent.min_link_travel_time[temp_link_no] += g_decelerating_time_slow_train;
									agent.max_link_travel_time[temp_link_no] += g_decelerating_time_slow_train;
								}
							}

						}
						else if (l == agent.set_of_allowed_links_LR.size() - 2)
						{
							//deceleration time
							if (agent.speed_grade == 1) //'G', fast train
							{
								agent.min_link_travel_time[temp_link_no] += g_decelerating_time_fast_train;
								agent.max_link_travel_time[temp_link_no] += g_decelerating_time_fast_train;
							}
							else if (agent.speed_grade == 2) //'D', slow train
							{
								agent.min_link_travel_time[temp_link_no] += g_decelerating_time_slow_train;
								agent.max_link_travel_time[temp_link_no] += g_decelerating_time_slow_train;
							}

							//acceleration time
							former_link_no = agent.set_of_allowed_links_LR[l - 1];

							if (g_link_vector[former_link_no].link_type == siding_track_type)
							{
								if (agent.speed_grade == 1) //'G', fast train
								{
									agent.min_link_travel_time[temp_link_no] += g_accelerating_time_fast_train;
									agent.max_link_travel_time[temp_link_no] += g_accelerating_time_fast_train;
								}
								else if (agent.speed_grade == 2) //'D', slow train
								{
									agent.min_link_travel_time[temp_link_no] += g_accelerating_time_slow_train;
									agent.max_link_travel_time[temp_link_no] += g_accelerating_time_slow_train;
								}
							}

						}
						else //intermediate section
						{
							former_link_no = agent.set_of_allowed_links_LR[l - 1];
							next_link_no = agent.set_of_allowed_links_LR[l + 1];

							//acceleration time
							if (g_link_vector[former_link_no].link_type == siding_track_type)
							{
								if (agent.speed_grade == 1) //'G', fast train
								{
									agent.min_link_travel_time[temp_link_no] += g_accelerating_time_fast_train;
									agent.max_link_travel_time[temp_link_no] += g_accelerating_time_fast_train;
								}
								else if (agent.speed_grade == 2) //'D', slow train
								{
									agent.min_link_travel_time[temp_link_no] += g_accelerating_time_slow_train;
									agent.max_link_travel_time[temp_link_no] += g_accelerating_time_slow_train;
								}
							}

							//deceleration time
							if (g_link_vector[next_link_no].link_type == siding_track_type)
							{
								if (agent.speed_grade == 1) //'G', fast train
								{
									agent.min_link_travel_time[temp_link_no] += g_decelerating_time_fast_train;
									agent.max_link_travel_time[temp_link_no] += g_decelerating_time_fast_train;
								}
								else if (agent.speed_grade == 2) //'D', slow train
								{
									agent.min_link_travel_time[temp_link_no] += g_decelerating_time_slow_train;
									agent.max_link_travel_time[temp_link_no] += g_decelerating_time_slow_train;
								}
							}
						}
					}
				}

				g_agent_vector.push_back(agent);
				g_number_of_agents++;
				if (g_number_of_agents % 1000 == 0)
					cout << "reading = " << g_number_of_agents / 1000 << " k agents..." << endl;
			}
		}

		parser_agent.CloseCSVFile();
	}

	cout << "number of agents = " << g_agent_vector.size() << endl;
}

//class for space time states
class CSTS_State
{
public:
	float m_speed;

	std::vector<int> m_outgoing_state_index_vector;
	std::vector<int> m_outgoing_state_transition_cost_vector;

};

std::vector<CSTS_State> g_STSStateVector;

void g_add_state_transition(int from_state, int to_state, float TransitionCost)
{
	g_STSStateVector[from_state].m_outgoing_state_index_vector.push_back(to_state);
	g_STSStateVector[from_state].m_outgoing_state_transition_cost_vector.push_back(TransitionCost);
}

class STSNetwork  // mainly for STS shortest path calculation
{
public:
	int m_threadNo;  // internal thread number 
	std::vector<int>  m_agent_vector; // assigned agents for computing 

	int m_number_of_nodes, m_number_of_time_intervals, m_number_of_states;

	int m_origin_node;
	int m_departure_time_beginning;
	int m_arrival_time_ending;

	float*** m_label_cost;
	int*** m_node_predecessor;
	int*** m_time_predecessor;
	int*** m_state_predecessor;

	STSNetwork()
	{
		m_origin_node = -1;
		m_label_cost = NULL;
		m_node_predecessor = NULL;
		m_time_predecessor = NULL;
		m_state_predecessor = NULL;
	}

	void AllocateSTSMemory(int number_of_nodes, int number_of_time_intervals, int number_of_states)
	{
		m_number_of_nodes = number_of_nodes;
		m_number_of_time_intervals = number_of_time_intervals;
		m_number_of_states = number_of_states;

		m_label_cost = Allocate3DDynamicArray<float>(m_number_of_nodes, m_number_of_time_intervals, m_number_of_states);
		m_node_predecessor = Allocate3DDynamicArray<int>(m_number_of_nodes, m_number_of_time_intervals, m_number_of_states);
		m_time_predecessor = Allocate3DDynamicArray<int>(m_number_of_nodes, m_number_of_time_intervals, m_number_of_states);
		m_state_predecessor = Allocate3DDynamicArray<int>(m_number_of_nodes, m_number_of_time_intervals, m_number_of_states);
	}

	~STSNetwork()
	{
		Deallocate3DDynamicArray<float>(m_label_cost, m_number_of_nodes, m_number_of_time_intervals);
		Deallocate3DDynamicArray<int>(m_node_predecessor, m_number_of_nodes, m_number_of_time_intervals);
		Deallocate3DDynamicArray<int>(m_time_predecessor, m_number_of_nodes, m_number_of_time_intervals);
		Deallocate3DDynamicArray<int>(m_state_predecessor, m_number_of_nodes, m_number_of_time_intervals);
	}

	//parallel computing version
	float optimal_STS_dynamic_programming(int departure_time_beginning, int arrival_time_ending, int deduce_flag)
	{
		if (m_origin_node < 0)
			return -1;

		float total_cost = _MAX_LABEL_COST;

		if (g_node_vector[m_origin_node].m_outgoing_node_vector.size() == 0)
		{
			return _MAX_LABEL_COST;
		}

		if (arrival_time_ending > m_number_of_time_intervals - 1)
		{
			return _MAX_LABEL_COST;
		}

		// step 1: Initialization for all nodes
		for (int i = 0; i < m_number_of_nodes; i++)
		{
			// to do: only update node label on the agent path
			for (int t = departure_time_beginning; t <= arrival_time_ending; t++)
			{
				for (int w = 0; w < m_number_of_states; w++)
				{
					m_label_cost[i][t][w] = _MAX_LABEL_COST;
					m_node_predecessor[i][t][w] = -1;  // pointer to previous NODE INDEX from the current label at current node and time
					m_time_predecessor[i][t][w] = -1;  // pointer to previous TIME INDEX from the current label at current node and time
					m_state_predecessor[i][t][w] = -1;
				}
			}
		}

		//step 2: Initialization for origin node at the preferred departure time, at departure time
		CAgent* p_agent = &(g_agent_vector[m_agent_vector[0]]); //the first agent is the current agent

		int w0 = 0;  // start from empty

		for (int t = departure_time_beginning; t <= min(departure_time_beginning + p_agent->departure_time_window, g_number_of_simulation_intervals - 1); t++)  //first loop: time
		{
			m_label_cost[m_origin_node][t][w0] = 0;
		}

		if (g_dp_algorithm_debug_flag == 1)
		{
			fprintf(g_pFileDebugLog, "****************starting of the DP****************\n");
			fprintf(g_pFileDebugLog, "agent_id = %d\n", p_agent->agent_id);
		}

		for (int t = departure_time_beginning; t <= arrival_time_ending; t++)  //first loop: time
		{
			if (g_dp_algorithm_debug_flag == 1)
			{
				fprintf(g_pFileDebugLog, "t = %d\n", t);
			}

			for (int n = 0; n < g_node_vector.size(); n++)
			{
				int temp_node_no = g_node_vector[n].node_seq_no;

				for (int link = 0; link < g_node_vector[temp_node_no].m_outgoing_node_vector.size(); link++)
				{
					int link_no = g_node_vector[temp_node_no].m_outgoing_node_vector[link].link_seq_no;

					if (p_agent->m_set_of_allowed_links_flag_LR[link_no] == 0 || g_link_vector[link_no].link_type == dummy_track_type)
					{
						continue;
					}

					int from_node = g_node_vector[temp_node_no].m_outgoing_node_vector[link].from_node_seq_no;
					int to_node = g_node_vector[temp_node_no].m_outgoing_node_vector[link].to_node_seq_no;

					if (g_dp_algorithm_debug_flag == 1)
					{
						fprintf(g_pFileDebugLog, "link_id = %d, from_node = %d, to_node = %d\n", g_link_vector[link_no].link_seq_no, g_internal_node_seq_no_to_node_id_map[from_node],
							g_internal_node_seq_no_to_node_id_map[to_node]);
					}

					for (int travel_time = p_agent->min_link_travel_time[link_no]; travel_time <= p_agent->max_link_travel_time[link_no]; travel_time++)
					{
						for (int w1 = 0; w1 < m_number_of_states; w1++) // for each state
						{
							if (g_dp_algorithm_debug_flag == 1)
							{
								fprintf(g_pFileDebugLog, "w1 = %d\n", w1);
							}

							if (m_label_cost[from_node][t][w1] < _MAX_LABEL_COST - 1)  //for feasible time-space point only
							{
								for (int w2_index = 0; w2_index < g_STSStateVector[w1].m_outgoing_state_index_vector.size(); w2_index++)
								{
									int travel_cost = travel_time;

									int w2 = g_STSStateVector[w1].m_outgoing_state_index_vector[w2_index];

									int new_to_node_arrival_time = min(t + travel_time, g_number_of_simulation_intervals - 1); // plus the waiting time on the link

									float temporary_label_cost = 0;

									float sum_of_multipliers = 0;

									if (g_dp_algorithm_debug_flag == 1)
									{
										fprintf(g_pFileDebugLog, "w2 = %d, travel_time = %d\n", w2, travel_time);
									}

									bool isInfeasible = false;
									int nPosition = 0;
									int TA_left_time = t;
									int TD_right_time = new_to_node_arrival_time;
									int TA_left_time_headway = 0;
									int TD_right_time_headway = 0;
									int next_link_no = 0;
									int former_link_no = 0;

									vector <int>::iterator iElement = find(p_agent->set_of_allowed_links_LR.begin(), p_agent->set_of_allowed_links_LR.begin(), link_no);

									nPosition = distance(p_agent->set_of_allowed_links_LR.begin(), iElement);

									if (g_link_vector[link_no].link_type != siding_track_type && g_link_vector[link_no].link_type != dummy_track_type) // 11 is the station track
									{
										if (nPosition == 0) //first section
										{
											TA_left_time_headway = g_departure_headway_stop;
											next_link_no = p_agent->set_of_allowed_links_LR[nPosition + 1];

											if (g_link_vector[next_link_no].link_type == siding_track_type) //a arrival event
											{
												TD_right_time_headway = g_arrival_headway_stop;
											}
											else
											{
												TD_right_time_headway = g_arrival_headway_passing;
											}
										}
										else if (nPosition == (p_agent->set_of_allowed_links_LR.size() - 2)) //last section
										{
											TD_right_time_headway = g_arrival_headway_stop;

											former_link_no = p_agent->set_of_allowed_links_LR[nPosition - 1];

											if (g_link_vector[former_link_no].link_type == siding_track_type) //a departure event
											{
												TA_left_time_headway = g_departure_headway_stop;
											}
											else //a passing event
											{
												TA_left_time_headway = g_departure_headway_passing;
											}
										}
										else if (nPosition > 0 && nPosition < (p_agent->set_of_allowed_links_LR.size() - 2))//intermediate section
										{
											next_link_no = p_agent->set_of_allowed_links_LR[nPosition + 1];

											if (g_link_vector[next_link_no].link_type == siding_track_type) //a arrival event
											{
												TD_right_time_headway = g_arrival_headway_stop;
											}
											else
											{
												TD_right_time_headway = g_arrival_headway_passing;
											}

											former_link_no = p_agent->set_of_allowed_links_LR[nPosition - 1];

											if (g_link_vector[former_link_no].link_type == siding_track_type) //a departure event
											{
												TA_left_time_headway = g_departure_headway_stop;
											}
											else //a passing event
											{
												TA_left_time_headway = g_departure_headway_passing;
											}
										}
									}

									if (g_link_vector[link_no].link_type != dummy_track_type && g_link_vector[link_no].link_type != siding_track_type)
									{
										for (int time = TA_left_time; time <= min(TA_left_time + TA_left_time_headway - 1, g_number_of_simulation_intervals - 1); time++)
										{
											for (int fre = 0; fre < p_agent->frequency; fre++)
											{
												for (int h = 0; h <= H; h++)
												{
													int temp_time = min(time + fre * ceil(g_cycle_length / p_agent->frequency) + h * g_cycle_length, g_number_of_simulation_intervals - 1);

													if (g_link_vector[link_no].departure_state_dependent_travel_time_matrix[temp_time][w1] < 0)
													{
														isInfeasible = true;
														travel_cost = _MAX_LABEL_COST;
														break;
													}

													if (g_link_vector[link_no].same_link_id != 1000)
													{
														int same_link_no = g_internal_link_no_map[g_link_vector[link_no].same_link_id];

														if (g_link_vector[same_link_no].departure_state_dependent_travel_time_matrix[temp_time][w1] < 0)
														{
															isInfeasible = true;
															travel_cost = _MAX_LABEL_COST;
															break;
														}
													}
												}
											}

											if (isInfeasible == true)
											{
												break;
											}
										}

										for (int time = TD_right_time; time <= min(TD_right_time + TD_right_time_headway - 1, g_number_of_simulation_intervals - 1); time++)
										{
											for (int fre = 0; fre < p_agent->frequency; fre++)
											{
												for (int h = 0; h <= H; h++)
												{
													int temp_time = min(time + fre * ceil(g_cycle_length / p_agent->frequency) + h * g_cycle_length, g_number_of_simulation_intervals - 1);

													if (g_link_vector[link_no].arrival_state_dependent_travel_time_matrix[temp_time][w1] < 0)
													{
														isInfeasible = true;
														travel_cost = _MAX_LABEL_COST;
														break;
													}

													if (g_link_vector[link_no].same_link_id != 1000)
													{
														int same_link_no = g_internal_link_no_map[g_link_vector[link_no].same_link_id];

														if (g_link_vector[same_link_no].arrival_state_dependent_travel_time_matrix[temp_time][w1] < 0)
														{
															isInfeasible = true;
															travel_cost = _MAX_LABEL_COST;
															break;
														}
													}
												}
											}

											if (isInfeasible == true)
											{
												break;
											}
										}
									}

									if (deduce_flag == 0) //Lagrangian relaxation
									{
										for (int time = TA_left_time; time <= min(TA_left_time + TA_left_time_headway - 1, g_number_of_simulation_intervals - 1); time++)
										{
											for (int f = 0; f < p_agent->frequency; f++)
											{
												for (int h = 0; h <= H; h++)
												{
													int time_temp = min(time + f * ceil(g_cycle_length / p_agent->frequency) + h * g_cycle_length, g_number_of_simulation_intervals - 1);
													sum_of_multipliers += g_link_vector[link_no].state_time_dependent_departure_LR_multiplier_matrix[time_temp][w1];
												}
											}
										}

										for (int time = TD_right_time; time <= min(TD_right_time + TD_right_time_headway - 1, g_number_of_simulation_intervals - 1); time++)
										{
											for (int f = 0; f < p_agent->frequency; f++)
											{
												for (int h = 0; h <= H; h++)
												{
													int time_temp = min(time + f * ceil(g_cycle_length / p_agent->frequency) + h * g_cycle_length, g_number_of_simulation_intervals - 1);
													sum_of_multipliers += g_link_vector[link_no].state_time_dependent_arrival_LR_multiplier_matrix[time_temp][w1];
												}
											}
										}
									}
									else if (deduce_flag == 2) // ADMM
									{
										for (int time = TA_left_time; time <= min(TA_left_time + TA_left_time_headway - 1, g_number_of_simulation_intervals - 1); time++)
										{
											for (int f = 0; f < p_agent->frequency; f++)
											{
												for (int h = 0; h <= H; h++)
												{
													int time_temp = min(time + f * ceil(g_cycle_length / p_agent->frequency) + h * g_cycle_length, g_number_of_simulation_intervals - 1);
													sum_of_multipliers += g_link_vector[link_no].state_dependent_time_dependent_departure_ADMM_multiplier_matrix[time_temp][w1];
												}
											}
										}

										for (int time = TD_right_time; time <= min(TD_right_time + TD_right_time_headway - 1, g_number_of_simulation_intervals - 1); time++)
										{
											for (int f = 0; f < p_agent->frequency; f++)
											{
												for (int h = 0; h <= H; h++)
												{
													int time_temp = min(time + f * ceil(g_cycle_length / p_agent->frequency) + h * g_cycle_length, g_number_of_simulation_intervals - 1);
													sum_of_multipliers += g_link_vector[link_no].state_dependent_time_dependent_arrival_ADMM_multiplier_matrix[time_temp][w1];
												}
											}
										}
									}

									temporary_label_cost = m_label_cost[from_node][t][w1] + travel_cost * p_agent->frequency * 2 + sum_of_multipliers;

									if (temporary_label_cost < m_label_cost[to_node][new_to_node_arrival_time][w2])
									{
										// update cost label and node/time predecessor
										m_label_cost[to_node][new_to_node_arrival_time][w2] = temporary_label_cost;
										m_node_predecessor[to_node][new_to_node_arrival_time][w2] = from_node;
										m_time_predecessor[to_node][new_to_node_arrival_time][w2] = t;
										m_state_predecessor[to_node][new_to_node_arrival_time][w2] = w1;
									}

								}
							}  // feasible vertex label cost
						}  // for all states
					} //for each waiting time
				} //for each outgoing link
			}
		} // for all time t

		if (g_dp_algorithm_debug_flag == 1)
		{
			fprintf(g_pFileDebugLog, "****************End of the DP****************\n");
		}

		return total_cost;
	}

	float find_STS_path_for_agents_assigned_for_this_thread(int number_of_threads, int assignment_iteration_no, int deduce_flag)
	{
		int reversed_path_node_sequence[_MAX_NUMBER_OF_TIME_INTERVALS] = { -1 };
		int reversed_path_time_sequence[_MAX_NUMBER_OF_TIME_INTERVALS] = { -1 };
		int reversed_path_state_sequence[_MAX_NUMBER_OF_TIME_INTERVALS] = { -1 };
		float reversed_path_cost_sequence[_MAX_NUMBER_OF_TIME_INTERVALS] = { -1 };

		int path_node_sequence[_MAX_NUMBER_OF_TIME_INTERVALS] = { -1 };
		int path_link_sequence[_MAX_NUMBER_OF_TIME_INTERVALS] = { -1 };
		int path_time_sequence[_MAX_NUMBER_OF_TIME_INTERVALS] = { -1 };
		int path_state_sequence[_MAX_NUMBER_OF_TIME_INTERVALS] = { -1 };
		float path_cost_sequence[_MAX_NUMBER_OF_TIME_INTERVALS] = { -1 };

		// perform one to all STS shortest path
		int return_value = optimal_STS_dynamic_programming(m_departure_time_beginning, m_arrival_time_ending, deduce_flag);

		float total_cost = _MAX_LABEL_COST;

		for (int i = 0; i < m_agent_vector.size(); i++)
		{
			CAgent* p_agent = &(g_agent_vector[m_agent_vector[i]]);

			p_agent->path_link_seq_no_vector_upper_bound.clear();  // reset;
			p_agent->path_timestamp_vector_upper_bound.clear();
			p_agent->path_node_id_vector_upper_bound.clear();
			p_agent->path_new_link_id_vector_upper_bound.clear();
			p_agent->path_link_TA_vector_upper_bound.clear();
			p_agent->path_link_TD_vector_upper_bound.clear();
			p_agent->path_new_node_id_vector_upper_bound.clear();
			p_agent->path_new_timestamp_vector_upper_bound.clear();

			if (return_value == -1)
			{
				fprintf(g_pFileDebugLog, "agent %d with can not find destination node,\n", i);
				continue;
			}

			int current_node_seq_no;
			int current_link_seq_no;

			//back trace from the destination node to find the shortest path from shortest path tree 
			int destination_node_seq_no = g_internal_node_seq_no_map[p_agent->destination_node_id];

			int min_cost_time_index = m_arrival_time_ending;
			int w = 0;
			total_cost = m_label_cost[destination_node_seq_no][min_cost_time_index][w];

			for (int t = m_arrival_time_ending; t >m_departure_time_beginning; t--)
			{
				if (m_label_cost[destination_node_seq_no][t][w] < total_cost)
				{
					min_cost_time_index = t;
					total_cost = m_label_cost[destination_node_seq_no][t][w];
				}
			}

			//backtrack to the origin (based on node and time predecessors)
			int	node_size = 0;
			reversed_path_node_sequence[node_size] = destination_node_seq_no; //record the first node backward, destination node
			reversed_path_time_sequence[node_size] = min_cost_time_index;
			reversed_path_state_sequence[node_size] = w;
			reversed_path_cost_sequence[node_size] = m_label_cost[destination_node_seq_no][min_cost_time_index][w];

			node_size++;

			int pred_node = m_node_predecessor[destination_node_seq_no][min_cost_time_index][w];
			int pred_time = m_time_predecessor[destination_node_seq_no][min_cost_time_index][w];
			int pred_state = m_state_predecessor[destination_node_seq_no][min_cost_time_index][w];

			while (pred_node != -1 && node_size < _MAX_NUMBER_OF_TIME_INTERVALS) // scan backward in the predecessor array of the shortest path calculation results
			{
				reversed_path_node_sequence[node_size] = pred_node;
				reversed_path_time_sequence[node_size] = pred_time;
				reversed_path_state_sequence[node_size] = pred_state;
				reversed_path_cost_sequence[node_size] = m_label_cost[pred_node][pred_time][pred_state];

				node_size++;

				//record current values of node and time predecessors, and update PredNode and PredTime
				int pred_node_record = pred_node;
				int pred_time_record = pred_time;
				int pred_state_record = pred_state;

				pred_node = m_node_predecessor[pred_node_record][pred_time_record][pred_state_record];
				pred_time = m_time_predecessor[pred_node_record][pred_time_record][pred_state_record];
				pred_state = m_state_predecessor[pred_node_record][pred_time_record][pred_state_record];
			}

			//reverse the node sequence 
			for (int n = 0; n < node_size; n++)
			{
				path_node_sequence[n] = reversed_path_node_sequence[node_size - n - 1];
				path_time_sequence[n] = reversed_path_time_sequence[node_size - n - 1];
				path_state_sequence[n] = reversed_path_state_sequence[node_size - n - 1];
				path_cost_sequence[n] = reversed_path_cost_sequence[node_size - n - 1];
			}

			for (int temp_i = 0; temp_i < node_size; temp_i++)  // for each node 
			{
				p_agent->path_node_id_vector_upper_bound.push_back(g_internal_node_seq_no_to_node_id_map[path_node_sequence[temp_i]]);
				p_agent->path_timestamp_vector_upper_bound.push_back(path_time_sequence[temp_i]);
			}

			for (int temp_i = 0; temp_i < node_size - 1; temp_i++)  // for each link
			{
				int link_no = g_GetLinkSeqNo(g_internal_node_seq_no_map[p_agent->path_node_id_vector_upper_bound[temp_i]], g_internal_node_seq_no_map[p_agent->path_node_id_vector_upper_bound[temp_i + 1]]);
				path_link_sequence[temp_i] = link_no;

				if (link_no == -1)
				{
					continue;
				}

				p_agent->path_link_seq_no_vector_upper_bound.push_back(link_no);
			}

			float travel_time_return_value = path_time_sequence[node_size - 1] - path_time_sequence[0];

			int path_number_of_nodes = node_size;
		}

		for (int a = 0; a < m_agent_vector.size(); a++)
		{
			CAgent* p_agent = &(g_agent_vector[m_agent_vector[a]]);

			int temp_node_id = p_agent->path_node_id_vector_upper_bound[0];
			int node_arrival_time = p_agent->path_timestamp_vector_upper_bound[0];
			int node_departure_time = p_agent->path_timestamp_vector_upper_bound[0];

			for (int n = 0; n < p_agent->path_node_id_vector_upper_bound.size() - 1; n++)
			{
				if (p_agent->path_node_id_vector_upper_bound[n + 1] == temp_node_id)
				{
					node_departure_time = p_agent->path_timestamp_vector_upper_bound[n + 1];
				}
				else
				{
					p_agent->path_new_node_id_vector_upper_bound.push_back(temp_node_id);
					p_agent->path_new_timestamp_vector_upper_bound.push_back(node_arrival_time);
					p_agent->path_new_timestamp_vector_upper_bound.push_back(node_departure_time);

					temp_node_id = p_agent->path_node_id_vector_upper_bound[n + 1];
					node_arrival_time = p_agent->path_timestamp_vector_upper_bound[n + 1];
					node_departure_time = p_agent->path_timestamp_vector_upper_bound[n + 1];
				}
			}

			p_agent->path_new_node_id_vector_upper_bound.push_back(temp_node_id);
			p_agent->path_new_timestamp_vector_upper_bound.push_back(node_arrival_time);
			p_agent->path_new_timestamp_vector_upper_bound.push_back(node_departure_time);
		}

		//get the new path link vector, TA and TD vector
		for (int a = 0; a < m_agent_vector.size(); a++)
		{
			CAgent* p_agent = &(g_agent_vector[m_agent_vector[a]]);

			int temp_number = 1;

			for (int n = 0; n < p_agent->path_new_node_id_vector_upper_bound.size() - 1; n++)
			{
				p_agent->path_link_TA_vector_upper_bound.push_back(p_agent->path_new_timestamp_vector_upper_bound[temp_number]);
				temp_number += 2;
			}

			temp_number = 3;

			for (int n = 1; n < p_agent->path_new_node_id_vector_upper_bound.size(); n++)
			{
				p_agent->path_link_TD_vector_upper_bound.push_back(p_agent->path_new_timestamp_vector_upper_bound[temp_number]);
				temp_number += 2;
			}

			for (int n = 0; n < p_agent->path_new_node_id_vector_upper_bound.size() - 1; n++)
			{
				int from_node_no = g_internal_node_seq_no_map[p_agent->path_new_node_id_vector_upper_bound[n]];
				int to_node_no = g_internal_node_seq_no_map[p_agent->path_new_node_id_vector_upper_bound[n + 1]];
				int temp_link_no = g_GetLinkSeqNo(from_node_no, to_node_no);
				p_agent->path_new_link_id_vector_upper_bound.push_back(temp_link_no);
			}
		}

		//scan the shortest path to compute the time_dependent link volume
		for (int a = 0; a < m_agent_vector.size(); a++)
		{
			CAgent* p_agent = &(g_agent_vector[m_agent_vector[a]]);
			int w1 = 0;

			for (int l = 0; l < p_agent->path_new_link_id_vector_upper_bound.size(); l++)  // for each link in the path of this agent
			{
				int link_seq_no = p_agent->path_new_link_id_vector_upper_bound[l];

				//mark the link travel time in the same direction
				int TA = p_agent->path_link_TA_vector_upper_bound[l];
				int TD = p_agent->path_link_TD_vector_upper_bound[l];
				int TA_left_time = TA;
				int TD_right_time = TD;
				int TA_left_time_headway = 0;
				int TD_right_time_headway = 0;
				int next_link_no = 0;
				int former_link_no = 0;

				if (g_link_vector[link_seq_no].link_type != siding_track_type && g_link_vector[link_seq_no].link_type != dummy_track_type)
				{
					if (l == 0) //first section
					{
						TA_left_time_headway = g_departure_headway_stop;
						next_link_no = p_agent->path_new_link_id_vector_upper_bound[l + 1];

						if (g_link_vector[next_link_no].link_type == siding_track_type) //a arrival event
						{
							TD_right_time_headway = g_arrival_headway_stop;
						}
						else
						{
							TD_right_time_headway = g_arrival_headway_passing;
						}
					}
					else if (l == p_agent->path_new_link_id_vector_upper_bound.size() - 1) //last section
					{
						TD_right_time_headway = g_arrival_headway_stop;

						former_link_no = p_agent->path_new_link_id_vector_upper_bound[l - 1];

						if (g_link_vector[former_link_no].link_type == siding_track_type) //a departure event
						{
							TA_left_time_headway = g_departure_headway_stop;
						}
						else //a passing event
						{
							TA_left_time_headway = g_departure_headway_passing;
						}
					}
					else //intermediate section
					{
						next_link_no = p_agent->path_new_link_id_vector_upper_bound[l + 1];

						if (g_link_vector[next_link_no].link_type == siding_track_type) //a arrival event
						{
							TD_right_time_headway = g_arrival_headway_stop;
						}
						else
						{
							TD_right_time_headway = g_arrival_headway_passing;
						}

						former_link_no = p_agent->path_new_link_id_vector_upper_bound[l - 1];

						if (g_link_vector[former_link_no].link_type == siding_track_type) //a departure event
						{
							TA_left_time_headway = g_departure_headway_stop;
						}
						else //a passing event
						{
							TA_left_time_headway = g_departure_headway_passing;
						}
					}

					int time = 0;
					int time_next_cycle = 0;
					int same_link_no = 0;

					for (int temp_time = TA_left_time; temp_time <= min(TA_left_time + TA_left_time_headway - 1, g_number_of_simulation_intervals - 1); temp_time++)
					{
						for (int fre = 0; fre < p_agent->frequency; fre++)
						{
							for (int h = 0; h <= H; h++)
							{
								time = min(temp_time + fre * ceil(g_cycle_length / p_agent->frequency) + h * g_cycle_length, g_number_of_simulation_intervals - 1);
								g_link_vector[link_seq_no].departure_state_dependent_travel_time_matrix[time][w1] = g_link_vector[link_seq_no].departure_state_dependent_travel_time_matrix[time][w1] - 1;  // minus -1

								if (g_link_vector[link_seq_no].same_link_id != 1000)
								{
									same_link_no = g_internal_link_no_map[g_link_vector[link_seq_no].same_link_id];
									g_link_vector[same_link_no].departure_state_dependent_travel_time_matrix[time][w1] = g_link_vector[same_link_no].departure_state_dependent_travel_time_matrix[time][w1] - 1;  // minus -1
								}
							}
						}
					}

					for (int temp_time = TD_right_time; temp_time <= min(TD_right_time + TD_right_time_headway - 1, g_number_of_simulation_intervals - 1); temp_time++)
					{
						for (int fre = 0; fre < p_agent->frequency; fre++)
						{
							for (int h = 0; h <= H; h++)
							{
								time = min(temp_time + fre * ceil(g_cycle_length / p_agent->frequency) + h * g_cycle_length, g_number_of_simulation_intervals - 1);
								g_link_vector[link_seq_no].arrival_state_dependent_travel_time_matrix[time][w1] = g_link_vector[link_seq_no].arrival_state_dependent_travel_time_matrix[time][w1] - 1;  // minus -1

								if (g_link_vector[link_seq_no].same_link_id != 1000)
								{
									same_link_no = g_internal_link_no_map[g_link_vector[link_seq_no].same_link_id];
									g_link_vector[same_link_no].arrival_state_dependent_travel_time_matrix[time][w1] = g_link_vector[same_link_no].arrival_state_dependent_travel_time_matrix[time][w1] - 1;  // minus -1
								}
							}
						}
					}

				}

			}
		}

		return total_cost;
	}

	//state-time-space dynamic programming for Lagrangian relaxation
	float optimal_STS_dynamic_programming_LR(int departure_time_beginning, int arrival_time_ending, int LR_iteration)
	{
		if (m_origin_node < 0)
			return -1;

		float total_cost = _MAX_LABEL_COST;

		if (g_node_vector[m_origin_node].m_outgoing_node_vector.size() == 0)
		{
			return _MAX_LABEL_COST;
		}

		if (arrival_time_ending > m_number_of_time_intervals - 1)
		{
			return _MAX_LABEL_COST;
		}

		// step 1: Initialization for all nodes
		for (int i = 0; i < m_number_of_nodes; i++)
		{
			// to do: only update node label on the agent path
			for (int t = departure_time_beginning; t < g_number_of_simulation_intervals; t++)
			{
				for (int w = 0; w < m_number_of_states; w++)
				{
					m_label_cost[i][t][w] = _MAX_LABEL_COST;
					m_node_predecessor[i][t][w] = -1;  // pointer to previous NODE INDEX from the current label at current node and time
					m_time_predecessor[i][t][w] = -1;  // pointer to previous TIME INDEX from the current label at current node and time
					m_state_predecessor[i][t][w] = -1;
				}
			}
		}

		//Initialization for origin node at the preferred departure time, at departure time

		CAgent* p_agent = &(g_agent_vector[m_agent_vector[0]]); //the first agent is the current agent

		int w0 = 0;  // start from empty

		for (int t = departure_time_beginning; t <= min(departure_time_beginning + p_agent->departure_time_window, g_number_of_simulation_intervals - 1); t++)  //first loop: [0, T)
		{
			m_label_cost[m_origin_node][t][w0] = 0;
		}

		if (g_dp_algorithm_debug_flag == 1)
		{
			fprintf(g_pFileDebugLog, "****************starting of the DP for LR****************\n");
			fprintf(g_pFileDebugLog, "agent_id = %d\n", p_agent->agent_id);
		}

		for (int t = departure_time_beginning; t <= arrival_time_ending; t++)  //first loop: time
		{
			if (g_dp_algorithm_debug_flag == 1)
			{
				fprintf(g_pFileDebugLog, "t = %d\n", t);
			}

			for (int n = 0; n < g_node_vector.size(); n++)
			{
				int temp_node_no = g_node_vector[n].node_seq_no;

				for (int link = 0; link < g_node_vector[temp_node_no].m_outgoing_node_vector.size(); link++)
				{
					int link_no = g_node_vector[temp_node_no].m_outgoing_node_vector[link].link_seq_no;

					//is the link exist in the set of allowed links 
					//change to vector 0 not allowed, 1 allowed
					if (p_agent->m_set_of_allowed_links_flag_LR[link_no] == 0)
					{
						continue;
					}

					int from_node = g_node_vector[temp_node_no].m_outgoing_node_vector[link].from_node_seq_no;
					int to_node = g_node_vector[temp_node_no].m_outgoing_node_vector[link].to_node_seq_no;

					if (g_dp_algorithm_debug_flag == 1)
					{
						fprintf(g_pFileDebugLog, "link_id = %d, from_node = %d, to_node = %d\n", g_link_vector[link_no].link_seq_no, g_internal_node_seq_no_to_node_id_map[from_node],
							g_internal_node_seq_no_to_node_id_map[to_node]);
					}

					for (int travel_time = p_agent->min_link_travel_time[link_no]; travel_time <= p_agent->max_link_travel_time[link_no]; travel_time++)
					{
						for (int w1 = 0; w1 < m_number_of_states; w1++) // for each state
						{
							if (g_dp_algorithm_debug_flag == 1)
							{
								fprintf(g_pFileDebugLog, "w1 = %d\n", w1);
							}

							if (m_label_cost[from_node][t][w1] < _MAX_LABEL_COST - 1)  //for feasible time-space point only
							{
								for (int w2_index = 0; w2_index < g_STSStateVector[w1].m_outgoing_state_index_vector.size(); w2_index++)
								{
									int travel_cost = travel_time;

									int w2 = g_STSStateVector[w1].m_outgoing_state_index_vector[w2_index];

									int new_to_node_arrival_time = min(t + travel_time, g_number_of_simulation_intervals - 1); // plus the waiting time on the link

									float temporary_label_cost = 0;
									float sum_of_multipliers = 0;

									int nPosition = 0;
									int TA_left_time = t;
									int TD_right_time = new_to_node_arrival_time;
									int TA_left_time_headway = 0;
									int TD_right_time_headway = 0;
									int next_link_no = 0;
									int former_link_no = 0;

									vector <int>::iterator iElement = find(p_agent->set_of_allowed_links_LR.begin(), p_agent->set_of_allowed_links_LR.begin(), link_no);

									nPosition = distance(p_agent->set_of_allowed_links_LR.begin(), iElement);

									if (g_link_vector[link_no].link_type != siding_track_type && g_link_vector[link_no].link_type != dummy_track_type) // 11 is the station track
									{
										if (nPosition == 0) //first section
										{
											TA_left_time_headway = g_departure_headway_stop;
											next_link_no = p_agent->set_of_allowed_links_LR[nPosition + 1];

											if (g_link_vector[next_link_no].link_type == siding_track_type) //a arrival event
											{
												TD_right_time_headway = g_arrival_headway_stop;
											}
											else
											{
												TD_right_time_headway = g_arrival_headway_passing;
											}
										}
										else if (nPosition == (p_agent->set_of_allowed_links_LR.size() - 2)) //last section
										{
											TD_right_time_headway = g_arrival_headway_stop;

											former_link_no = p_agent->set_of_allowed_links_LR[nPosition - 1];

											if (g_link_vector[former_link_no].link_type == siding_track_type) //a departure event
											{
												TA_left_time_headway = g_departure_headway_stop;
											}
											else //a passing event
											{
												TA_left_time_headway = g_departure_headway_passing;
											}
										}
										else //intermediate section
										{
											next_link_no = p_agent->set_of_allowed_links_LR[nPosition + 1];

											if (g_link_vector[next_link_no].link_type == siding_track_type) //a arrival event
											{
												TD_right_time_headway = g_arrival_headway_stop;
											}
											else
											{
												TD_right_time_headway = g_arrival_headway_passing;
											}

											former_link_no = p_agent->set_of_allowed_links_LR[nPosition - 1];

											if (g_link_vector[former_link_no].link_type == siding_track_type) //a departure event
											{
												TA_left_time_headway = g_departure_headway_stop;
											}
											else //a passing event
											{
												TA_left_time_headway = g_departure_headway_passing;
											}
										}

										for (int time = TA_left_time; time <= min(TA_left_time + TA_left_time_headway - 1, g_number_of_simulation_intervals - 1); time++)
										{
											for (int f = 0; f < p_agent->frequency; f++)
											{
												for (int h = 0; h <= H; h++)
												{
													int time_temp = min(time + f * ceil(g_cycle_length / p_agent->frequency) + h * g_cycle_length, g_number_of_simulation_intervals - 1);
													sum_of_multipliers += g_link_vector[link_no].state_time_dependent_departure_LR_multiplier_matrix[time_temp][w1];
												}
											}
										}

										for (int time = TD_right_time; time <= min(TD_right_time + TD_right_time_headway - 1, g_number_of_simulation_intervals - 1); time++)
										{
											for (int f = 0; f < p_agent->frequency; f++)
											{
												for (int h = 0; h <= H; h++)
												{
													int time_temp = min(time + f * ceil(g_cycle_length / p_agent->frequency) + h * g_cycle_length, g_number_of_simulation_intervals - 1);
													sum_of_multipliers += g_link_vector[link_no].state_time_dependent_arrival_LR_multiplier_matrix[time_temp][w1];
												}
											}
										}
									}

									temporary_label_cost = m_label_cost[from_node][t][w1] + travel_cost * p_agent->frequency * 2 + sum_of_multipliers;

									if (temporary_label_cost < m_label_cost[to_node][new_to_node_arrival_time][w2]) // we only compare cost at the downstream node ToID at the new arrival time t
									{
										// update cost label and node/time predecessor
										m_label_cost[to_node][new_to_node_arrival_time][w2] = temporary_label_cost;
										m_node_predecessor[to_node][new_to_node_arrival_time][w2] = from_node;  // pointer to previous NODE INDEX from the current label at current node and time
										m_time_predecessor[to_node][new_to_node_arrival_time][w2] = t;  // pointer to previous TIME INDEX from the current label at current node and time
										m_state_predecessor[to_node][new_to_node_arrival_time][w2] = w1;
									}
								}
							}  // feasible vertex label cost
						}  // for all states
					}
				} //for each outgoing link
			}
		} // for all time t

		if (g_dp_algorithm_debug_flag == 1)
		{
			fprintf(g_pFileDebugLog, "****************End of the DP for LR****************\n");
		}

		return total_cost;
	}

	float find_STS_path_for_agents_assigned_for_this_thread_LR(int number_of_threads, int LR_iteration)
	{
		int reversed_path_node_sequence[_MAX_NUMBER_OF_TIME_INTERVALS] = { -1 };
		int reversed_path_time_sequence[_MAX_NUMBER_OF_TIME_INTERVALS] = { -1 };
		int reversed_path_state_sequence[_MAX_NUMBER_OF_TIME_INTERVALS] = { -1 };
		float reversed_path_cost_sequence[_MAX_NUMBER_OF_TIME_INTERVALS] = { -1 };

		int path_node_sequence[_MAX_NUMBER_OF_TIME_INTERVALS] = { -1 };
		int path_link_sequence[_MAX_NUMBER_OF_TIME_INTERVALS] = { -1 };
		int path_time_sequence[_MAX_NUMBER_OF_TIME_INTERVALS] = { -1 };
		int path_state_sequence[_MAX_NUMBER_OF_TIME_INTERVALS] = { -1 };
		float path_cost_sequence[_MAX_NUMBER_OF_TIME_INTERVALS] = { -1 };

		// perform one to all STS shortest path
		int return_value = optimal_STS_dynamic_programming_LR(m_departure_time_beginning, m_arrival_time_ending, LR_iteration);

		float total_cost = _MAX_LABEL_COST;

		for (int i = 0; i < m_agent_vector.size(); i++)
		{
			CAgent* p_agent = &(g_agent_vector[m_agent_vector[i]]);

			p_agent->path_link_seq_no_vector.clear();  // reset;
			p_agent->path_timestamp_vector.clear();
			p_agent->path_node_id_vector.clear();

			if (return_value == -1)
			{
				fprintf(g_pFileDebugLog, "agent %d with can not find destination node,\n", i);
				continue;
			}

			int current_node_seq_no;
			int current_link_seq_no;

			//back trace from the destination node to find the shortest path from shortest path tree
			int destination_node_seq_no = g_internal_node_seq_no_map[p_agent->destination_node_id];

			int min_cost_time_index = m_arrival_time_ending;
			int w = 0;
			total_cost = m_label_cost[destination_node_seq_no][min_cost_time_index][w];

			for (int t = m_arrival_time_ending; t >m_departure_time_beginning; t--)
			{
				if (m_label_cost[destination_node_seq_no][t][w] <= total_cost)
				{
					min_cost_time_index = t;
					total_cost = m_label_cost[destination_node_seq_no][t][w];
				}
			}

			//backtrack to the origin (based on node and time predecessors)
			int	node_size = 0;
			reversed_path_node_sequence[node_size] = destination_node_seq_no; //record the first node backward, destination node
			reversed_path_time_sequence[node_size] = min_cost_time_index;
			reversed_path_state_sequence[node_size] = w;
			reversed_path_cost_sequence[node_size] = m_label_cost[destination_node_seq_no][min_cost_time_index][w];

			node_size++;

			int pred_node = m_node_predecessor[destination_node_seq_no][min_cost_time_index][w];
			int pred_time = m_time_predecessor[destination_node_seq_no][min_cost_time_index][w];
			int pred_state = m_state_predecessor[destination_node_seq_no][min_cost_time_index][w];

			while (pred_node != -1 && node_size < _MAX_NUMBER_OF_TIME_INTERVALS) // scan backward in the predecessor array of the shortest path calculation results
			{
				reversed_path_node_sequence[node_size] = pred_node;
				reversed_path_time_sequence[node_size] = pred_time;
				reversed_path_state_sequence[node_size] = pred_state;
				reversed_path_cost_sequence[node_size] = m_label_cost[pred_node][pred_time][pred_state];

				node_size++;

				//record current values of node and time predecessors, and update PredNode and PredTime
				int pred_node_record = pred_node;
				int pred_time_record = pred_time;
				int pred_state_record = pred_state;

				pred_node = m_node_predecessor[pred_node_record][pred_time_record][pred_state_record];
				pred_time = m_time_predecessor[pred_node_record][pred_time_record][pred_state_record];
				pred_state = m_state_predecessor[pred_node_record][pred_time_record][pred_state_record];
			}

			//reverse the node sequence 
			for (int n = 0; n < node_size; n++)
			{
				path_node_sequence[n] = reversed_path_node_sequence[node_size - n - 1];
				path_time_sequence[n] = reversed_path_time_sequence[node_size - n - 1];
				path_state_sequence[n] = reversed_path_state_sequence[node_size - n - 1];
				path_cost_sequence[n] = reversed_path_cost_sequence[node_size - n - 1];
			}

			for (int temp_i = 0; temp_i < node_size; temp_i++)  // for each node 
			{
				p_agent->path_node_id_vector.push_back(g_internal_node_seq_no_to_node_id_map[path_node_sequence[temp_i]]);
				p_agent->path_timestamp_vector.push_back(path_time_sequence[temp_i]);
			}

			for (int temp_i = 0; temp_i < node_size - 1; temp_i++)  // for each link, 
			{
				int link_no = g_GetLinkSeqNo(g_internal_node_seq_no_map[p_agent->path_node_id_vector[temp_i]], g_internal_node_seq_no_map[p_agent->path_node_id_vector[temp_i + 1]]);
				path_link_sequence[temp_i] = link_no;

				if (link_no == -1)
				{
					continue;
				}

				p_agent->path_link_seq_no_vector.push_back(link_no);
			}

			float travel_time_return_value = path_time_sequence[node_size - 1] - path_time_sequence[0];

			int path_number_of_nodes = node_size;
		}

		//get the arrival time and departure time for each node	
		for (int a = 0; a < m_agent_vector.size(); a++)
		{
			CAgent* p_agent = &(g_agent_vector[m_agent_vector[a]]);

			int temp_node_id = p_agent->path_node_id_vector[0];
			int node_arrival_time = p_agent->path_timestamp_vector[0];
			int node_departure_time = p_agent->path_timestamp_vector[0];

			for (int n = 0; n < p_agent->path_node_id_vector.size() - 1; n++)
			{
				if (p_agent->path_node_id_vector[n + 1] == temp_node_id)
				{
					node_departure_time = p_agent->path_timestamp_vector[n + 1];
				}
				else
				{
					p_agent->path_new_node_id_vector.push_back(temp_node_id);
					p_agent->path_new_timestamp_vector.push_back(node_arrival_time);
					p_agent->path_new_timestamp_vector.push_back(node_departure_time);

					temp_node_id = p_agent->path_node_id_vector[n + 1];
					node_arrival_time = p_agent->path_timestamp_vector[n + 1];
					node_departure_time = p_agent->path_timestamp_vector[n + 1];
				}
			}

			p_agent->path_new_node_id_vector.push_back(temp_node_id);
			p_agent->path_new_timestamp_vector.push_back(node_arrival_time);
			p_agent->path_new_timestamp_vector.push_back(node_departure_time);
		}

		//get the new path link vector, TA and TD vector
		for (int a = 0; a < m_agent_vector.size(); a++)
		{
			CAgent* p_agent = &(g_agent_vector[m_agent_vector[a]]);

			int temp_number = 1;

			for (int n = 0; n < p_agent->path_new_node_id_vector.size() - 1; n++)
			{
				p_agent->path_link_TA_vector.push_back(p_agent->path_new_timestamp_vector[temp_number]);
				temp_number += 2;
			}

			temp_number = 3;

			for (int n = 1; n < p_agent->path_new_node_id_vector.size(); n++)
			{
				p_agent->path_link_TD_vector.push_back(p_agent->path_new_timestamp_vector[temp_number]);
				temp_number += 2;
			}

			for (int n = 0; n < p_agent->path_new_node_id_vector.size() - 1; n++)
			{
				int from_node_no = g_internal_node_seq_no_map[p_agent->path_new_node_id_vector[n]];
				int to_node_no = g_internal_node_seq_no_map[p_agent->path_new_node_id_vector[n + 1]];
				int temp_link_no = g_GetLinkSeqNo(from_node_no, to_node_no);

				p_agent->path_new_link_id_vector.push_back(temp_link_no);
			}
		}

		CAgent* p_current_agent = &(g_agent_vector[m_agent_vector[0]]); //the first agent is the current agent

		int n = p_current_agent->path_link_TD_vector.size() - 1;

		if (LR_iteration == 0)
		{
			p_current_agent->free_flow_travel_time = p_current_agent->path_link_TD_vector[n] - p_current_agent->path_link_TA_vector[0];
			p_current_agent->travel_time_in_dual_solution = p_current_agent->path_link_TD_vector[n] - p_current_agent->path_link_TA_vector[0];
		}
		else
		{
			p_current_agent->travel_time_in_dual_solution = p_current_agent->path_link_TD_vector[n] - p_current_agent->path_link_TA_vector[0];
		}

		return total_cost;
	}

	//state-time-space dynamic programming for ADMM
	float optimal_STS_dynamic_programming_ADMM(int departure_time_beginning, int arrival_time_ending, int ADMM_iteration)
	{
		if (m_origin_node < 0)
			return -1;

		float total_cost = _MAX_LABEL_COST;

		if (g_node_vector[m_origin_node].m_outgoing_node_vector.size() == 0)
		{
			return _MAX_LABEL_COST;
		}

		if (arrival_time_ending > m_number_of_time_intervals - 1)
		{
			return _MAX_LABEL_COST;
		}

		// step 1: Initialization for all nodes
		for (int i = 0; i < m_number_of_nodes; i++)
		{
			// to do: only update node label on the agent path
			for (int t = departure_time_beginning; t < g_number_of_simulation_intervals; t++)
			{
				for (int w = 0; w < m_number_of_states; w++)
				{
					m_label_cost[i][t][w] = _MAX_LABEL_COST;
					m_node_predecessor[i][t][w] = -1;  // pointer to previous NODE INDEX from the current label at current node and time
					m_time_predecessor[i][t][w] = -1;  // pointer to previous TIME INDEX from the current label at current node and time
					m_state_predecessor[i][t][w] = -1;
				}
			}
		}

		CAgent* p_agent = &(g_agent_vector[m_agent_vector[0]]); //the first agent is the current agent

																//step 2: Initialization for origin node at the preferred departure time, at departure time

		int w0 = 0;  // start from empty

		for (int t = departure_time_beginning; t <= min(departure_time_beginning + p_agent->departure_time_window, g_number_of_simulation_intervals - 1); t++)  //first loop: [0, T)
		{
			m_label_cost[m_origin_node][t][w0] = 0;
		}

		//dynamic programming
		if (g_dp_algorithm_debug_flag == 1)
		{
			fprintf(g_pFileDebugLog, "****************starting of the DP for LR****************\n");
			fprintf(g_pFileDebugLog, "agent_id = %d\n", p_agent->agent_id);
		}

		for (int t = departure_time_beginning; t <= arrival_time_ending; t++)  //first loop: time
		{
			if (g_dp_algorithm_debug_flag == 1)
			{
				fprintf(g_pFileDebugLog, "t = %d\n", t);
			}

			for (int n = 0; n < g_node_vector.size(); n++)
			{
				int temp_node_no = g_node_vector[n].node_seq_no;

				for (int link = 0; link < g_node_vector[temp_node_no].m_outgoing_node_vector.size(); link++)
				{
					int link_no = g_node_vector[temp_node_no].m_outgoing_node_vector[link].link_seq_no;

					//is the link exist in the set of allowed links 
					//change to vector 0 not allowed, 1 allowed
					if (p_agent->m_set_of_allowed_links_flag_LR[link_no] == 0)
					{
						continue;
					}

					int from_node = g_node_vector[temp_node_no].m_outgoing_node_vector[link].from_node_seq_no;
					int to_node = g_node_vector[temp_node_no].m_outgoing_node_vector[link].to_node_seq_no;

					if (g_dp_algorithm_debug_flag == 1)
					{
						fprintf(g_pFileDebugLog, "link_id = %d, from_node = %d, to_node = %d\n", g_link_vector[link_no].link_seq_no, g_internal_node_seq_no_to_node_id_map[from_node],
							g_internal_node_seq_no_to_node_id_map[to_node]);
					}

					for (int travel_time = p_agent->min_link_travel_time[link_no]; travel_time <= p_agent->max_link_travel_time[link_no]; travel_time++)
					{
						for (int w1 = 0; w1 < m_number_of_states; w1++) // for each state
						{
							if (g_dp_algorithm_debug_flag == 1)
							{
								fprintf(g_pFileDebugLog, "w1 = %d\n", w1);
							}

							if (m_label_cost[from_node][t][w1] < _MAX_LABEL_COST - 1)  //for feasible time-space point only
							{
								for (int w2_index = 0; w2_index < g_STSStateVector[w1].m_outgoing_state_index_vector.size(); w2_index++)
								{
									int travel_cost = travel_time;

									int w2 = g_STSStateVector[w1].m_outgoing_state_index_vector[w2_index];

									int new_to_node_arrival_time = min(t + travel_time, g_number_of_simulation_intervals - 1); // plus the waiting time on the link

									float temporary_label_cost = 0;

									float sum_of_multipliers = 0;

									int nPosition = 0;
									int TA_left_time = t;
									int TD_right_time = new_to_node_arrival_time;
									int TA_left_time_headway = 0;
									int TD_right_time_headway = 0;
									int next_link_no = 0;
									int former_link_no = 0;

									vector <int>::iterator iElement = find(p_agent->set_of_allowed_links_LR.begin(), p_agent->set_of_allowed_links_LR.begin(), link_no);

									nPosition = distance(p_agent->set_of_allowed_links_LR.begin(), iElement);

									if (g_link_vector[link_no].link_type != siding_track_type && g_link_vector[link_no].link_type != dummy_track_type)
									{
										if (nPosition == 0) //first section
										{
											TA_left_time_headway = g_departure_headway_stop;
											next_link_no = p_agent->set_of_allowed_links_LR[nPosition + 1];

											if (g_link_vector[next_link_no].link_type == siding_track_type) //a arrival event
											{
												TD_right_time_headway = g_arrival_headway_stop;
											}
											else
											{
												TD_right_time_headway = g_arrival_headway_passing;
											}
										}
										else if (nPosition == (p_agent->set_of_allowed_links_LR.size() - 2)) //last section
										{
											TD_right_time_headway = g_arrival_headway_stop;

											former_link_no = p_agent->set_of_allowed_links_LR[nPosition - 1];

											if (g_link_vector[former_link_no].link_type == siding_track_type) //a departure event
											{
												TA_left_time_headway = g_departure_headway_stop;
											}
											else //a passing event
											{
												TA_left_time_headway = g_departure_headway_passing;
											}
										}
										else //intermediate section
										{
											next_link_no = p_agent->set_of_allowed_links_LR[nPosition + 1];

											if (g_link_vector[next_link_no].link_type == siding_track_type) //a arrival event
											{
												TD_right_time_headway = g_arrival_headway_stop;
											}
											else
											{
												TD_right_time_headway = g_arrival_headway_passing;
											}

											former_link_no = p_agent->set_of_allowed_links_LR[nPosition - 1];

											if (g_link_vector[former_link_no].link_type == siding_track_type) //a departure event
											{
												TA_left_time_headway = g_departure_headway_stop;
											}
											else //a passing event
											{
												TA_left_time_headway = g_departure_headway_passing;
											}
										}

										for (int time = TA_left_time; time <= min(TA_left_time + TA_left_time_headway - 1, g_number_of_simulation_intervals - 1); time++)
										{
											for (int f = 0; f < p_agent->frequency; f++)
											{
												for (int h = 0; h <= H; h++)
												{
													int time_temp = min(time + f * ceil(g_cycle_length / p_agent->frequency) + h * g_cycle_length, g_number_of_simulation_intervals - 1);
													sum_of_multipliers += g_link_vector[link_no].state_dependent_time_dependent_departure_ADMM_multiplier_matrix[time_temp][w1];
												}
											}
										}

										for (int time = TD_right_time; time <= min(TD_right_time + TD_right_time_headway - 1, g_number_of_simulation_intervals - 1); time++)
										{
											for (int f = 0; f < p_agent->frequency; f++)
											{
												for (int h = 0; h <= H; h++)
												{
													int time_temp = min(time + f * ceil(g_cycle_length / p_agent->frequency) + h * g_cycle_length, g_number_of_simulation_intervals - 1);
													sum_of_multipliers += g_link_vector[link_no].state_dependent_time_dependent_arrival_ADMM_multiplier_matrix[time_temp][w1];
												}
											}
										}
									}

									temporary_label_cost = m_label_cost[from_node][t][w1] + travel_cost * p_agent->frequency * 2 + sum_of_multipliers;

									if (temporary_label_cost < m_label_cost[to_node][new_to_node_arrival_time][w2]) // we only compare cost at the downstream node ToID at the new arrival time t
									{
										// update cost label and node/time predecessor
										m_label_cost[to_node][new_to_node_arrival_time][w2] = temporary_label_cost;
										m_node_predecessor[to_node][new_to_node_arrival_time][w2] = from_node;  // pointer to previous NODE INDEX from the current label at current node and time
										m_time_predecessor[to_node][new_to_node_arrival_time][w2] = t;  // pointer to previous TIME INDEX from the current label at current node and time
										m_state_predecessor[to_node][new_to_node_arrival_time][w2] = w1;
									}
								}
							}  // feasible vertex label cost
						}  // for all states
					}

				} //for each outgoing link
			}
		} // for all time t

		if (g_dp_algorithm_debug_flag == 1)
		{
			fprintf(g_pFileDebugLog, "****************End of the DP for ADMM****************\n");
		}

		return total_cost;
	}

	float find_STS_path_for_agents_assigned_for_this_thread_ADMM(int number_of_threads, int ADMM_iteration)
	{
		int reversed_path_node_sequence[_MAX_NUMBER_OF_TIME_INTERVALS] = { -1 };
		int reversed_path_time_sequence[_MAX_NUMBER_OF_TIME_INTERVALS] = { -1 };
		int reversed_path_state_sequence[_MAX_NUMBER_OF_TIME_INTERVALS] = { -1 };
		float reversed_path_cost_sequence[_MAX_NUMBER_OF_TIME_INTERVALS] = { -1 };

		int path_node_sequence[_MAX_NUMBER_OF_TIME_INTERVALS] = { -1 };
		int path_link_sequence[_MAX_NUMBER_OF_TIME_INTERVALS] = { -1 };
		int path_time_sequence[_MAX_NUMBER_OF_TIME_INTERVALS] = { -1 };
		int path_state_sequence[_MAX_NUMBER_OF_TIME_INTERVALS] = { -1 };
		float path_cost_sequence[_MAX_NUMBER_OF_TIME_INTERVALS] = { -1 };

		// perform one to all state-time-state shortest path
		int return_value = optimal_STS_dynamic_programming_ADMM(m_departure_time_beginning, m_arrival_time_ending, ADMM_iteration);

		float total_cost = _MAX_LABEL_COST;

		for (int i = 0; i < m_agent_vector.size(); i++)
		{
			CAgent* p_agent = &(g_agent_vector[m_agent_vector[i]]);

			p_agent->path_link_seq_no_vector.clear();  // reset;
			p_agent->path_timestamp_vector.clear();
			p_agent->path_node_id_vector.clear();
			p_agent->path_new_link_id_vector.clear();
			p_agent->path_link_TA_vector.clear();
			p_agent->path_link_TD_vector.clear();
			p_agent->path_new_node_id_vector.clear();
			p_agent->path_new_timestamp_vector.clear();
			p_agent->partial_schedule_node_vector.clear();
			p_agent->partial_schedule_node_TD.clear();
			p_agent->partial_schedule_node_TA.clear();

			if (return_value == -1)
			{
				fprintf(g_pFileDebugLog, "agent %d with can not find destination node,\n", i);
				continue;
			}

			int current_node_seq_no;
			int current_link_seq_no;

			//back trace from the destination node to find the shortest path from shortest path tree
			int destination_node_seq_no = g_internal_node_seq_no_map[p_agent->destination_node_id];

			int min_cost_time_index = m_arrival_time_ending;
			int w = 0;
			total_cost = m_label_cost[destination_node_seq_no][min_cost_time_index][w];

			for (int t = m_arrival_time_ending; t >m_departure_time_beginning; t--)
			{
				if (m_label_cost[destination_node_seq_no][t][w] <= total_cost)
				{
					min_cost_time_index = t;
					total_cost = m_label_cost[destination_node_seq_no][t][w];
				}
			}

			//backtrack to the origin (based on node and time predecessors)
			int	node_size = 0;
			reversed_path_node_sequence[node_size] = destination_node_seq_no; //record the first node backward, destination node
			reversed_path_time_sequence[node_size] = min_cost_time_index;
			reversed_path_state_sequence[node_size] = w;
			reversed_path_cost_sequence[node_size] = m_label_cost[destination_node_seq_no][min_cost_time_index][w];

			node_size++;

			int pred_node = m_node_predecessor[destination_node_seq_no][min_cost_time_index][w];
			int pred_time = m_time_predecessor[destination_node_seq_no][min_cost_time_index][w];
			int pred_state = m_state_predecessor[destination_node_seq_no][min_cost_time_index][w];

			while (pred_node != -1 && node_size < _MAX_NUMBER_OF_TIME_INTERVALS) // scan backward in the predecessor array of the shortest path calculation results
			{
				reversed_path_node_sequence[node_size] = pred_node;
				reversed_path_time_sequence[node_size] = pred_time;
				reversed_path_state_sequence[node_size] = pred_state;
				reversed_path_cost_sequence[node_size] = m_label_cost[pred_node][pred_time][pred_state];

				node_size++;

				//record current values of node and time predecessors, and update PredNode and PredTime

				int pred_node_record = pred_node;
				int pred_time_record = pred_time;
				int pred_state_record = pred_state;

				pred_node = m_node_predecessor[pred_node_record][pred_time_record][pred_state_record];
				pred_time = m_time_predecessor[pred_node_record][pred_time_record][pred_state_record];
				pred_state = m_state_predecessor[pred_node_record][pred_time_record][pred_state_record];
			}

			//reverse the node sequence 
			for (int n = 0; n < node_size; n++)
			{
				path_node_sequence[n] = reversed_path_node_sequence[node_size - n - 1];
				path_time_sequence[n] = reversed_path_time_sequence[node_size - n - 1];
				path_state_sequence[n] = reversed_path_state_sequence[node_size - n - 1];
				path_cost_sequence[n] = reversed_path_cost_sequence[node_size - n - 1];
			}

			for (int temp_i = 0; temp_i < node_size; temp_i++)  // for each node 
			{
				p_agent->path_node_id_vector.push_back(g_internal_node_seq_no_to_node_id_map[path_node_sequence[temp_i]]);
				p_agent->path_timestamp_vector.push_back(path_time_sequence[temp_i]);
			}

			for (int temp_i = 0; temp_i < node_size - 1; temp_i++)  // for each link, 
			{
				int link_no = g_GetLinkSeqNo(g_internal_node_seq_no_map[p_agent->path_node_id_vector[temp_i]], g_internal_node_seq_no_map[p_agent->path_node_id_vector[temp_i + 1]]);
				path_link_sequence[temp_i] = link_no;

				if (link_no == -1)
				{
					continue;
				}

				p_agent->path_link_seq_no_vector.push_back(link_no);
			}

			float travel_time_return_value = path_time_sequence[node_size - 1] - path_time_sequence[0];

			int path_number_of_nodes = node_size;
		}

		//get the arrival time and departure time for each node	
		for (int a = 0; a < m_agent_vector.size(); a++)
		{
			CAgent* p_agent = &(g_agent_vector[m_agent_vector[a]]);

			int temp_node_id = p_agent->path_node_id_vector[0];
			int node_arrival_time = p_agent->path_timestamp_vector[0];
			int node_departure_time = p_agent->path_timestamp_vector[0];

			for (int n = 0; n < p_agent->path_node_id_vector.size() - 1; n++)
			{
				if (p_agent->path_node_id_vector[n + 1] == temp_node_id)
				{
					node_departure_time = p_agent->path_timestamp_vector[n + 1];
				}
				else
				{
					p_agent->path_new_node_id_vector.push_back(temp_node_id);
					p_agent->path_new_timestamp_vector.push_back(node_arrival_time);
					p_agent->path_new_timestamp_vector.push_back(node_departure_time);

					temp_node_id = p_agent->path_node_id_vector[n + 1];
					node_arrival_time = p_agent->path_timestamp_vector[n + 1];
					node_departure_time = p_agent->path_timestamp_vector[n + 1];
				}
			}

			p_agent->path_new_node_id_vector.push_back(temp_node_id); //push back the last node
			p_agent->path_new_timestamp_vector.push_back(node_arrival_time);
			p_agent->path_new_timestamp_vector.push_back(node_departure_time);
		}

		//get the new path link vector, TA and TD vector
		for (int a = 0; a < m_agent_vector.size(); a++)
		{
			CAgent* p_agent = &(g_agent_vector[m_agent_vector[a]]);

			int temp_number = 1;

			for (int n = 0; n < p_agent->path_new_node_id_vector.size() - 1; n++)
			{
				p_agent->path_link_TA_vector.push_back(p_agent->path_new_timestamp_vector[temp_number]);
				temp_number += 2;
			}

			temp_number = 3;

			for (int n = 1; n < p_agent->path_new_node_id_vector.size(); n++)
			{
				p_agent->path_link_TD_vector.push_back(p_agent->path_new_timestamp_vector[temp_number]);
				temp_number += 2;
			}

			for (int n = 0; n < p_agent->path_new_node_id_vector.size() - 1; n++)
			{
				int from_node_no = g_internal_node_seq_no_map[p_agent->path_new_node_id_vector[n]];
				int to_node_no = g_internal_node_seq_no_map[p_agent->path_new_node_id_vector[n + 1]];
				int temp_link_no = g_GetLinkSeqNo(from_node_no, to_node_no);
				p_agent->path_new_link_id_vector.push_back(temp_link_no);
			}
		}

		CAgent* p_current_agent = &(g_agent_vector[m_agent_vector[0]]); //the first agent is the current agent

		int n = p_current_agent->path_link_TD_vector.size() - 1;

		if (ADMM_iteration == 0 || g_reoptimization_flag == true)
		{
			p_current_agent->free_flow_travel_time = p_current_agent->path_link_TD_vector[n] - p_current_agent->path_link_TA_vector[0];
			p_current_agent->travel_time_in_dual_solution = p_current_agent->path_link_TD_vector[n] - p_current_agent->path_link_TA_vector[0];
		}
		else
		{
			p_current_agent->travel_time_in_dual_solution = p_current_agent->path_link_TD_vector[n] - p_current_agent->path_link_TA_vector[0];
		}

		return total_cost;
	}
};

int g_number_of_CPU_threads()
{
	int number_of_threads = omp_get_max_threads();

	int max_number_of_threads = 8;

	if (number_of_threads > max_number_of_threads)
		number_of_threads = max_number_of_threads;

	return number_of_threads;
}

STSNetwork* pSTSNetwork = NULL;

void g_OutAgentCSVFile_FromSimulation_LR()
{
	FILE* g_pFileAgent_LR = NULL;

	int train_id = 0;
	int line_id = 0;
	int station_id = 0;
	std::string station_name;
	int x_time_in_min = 0;
	float y_station_location = 0;
	int g_draw_seq_no = 0;
	int temp_y_station_location = 0;

	g_pFileAgent_LR = fopen("output_agent_LR.csv", "w");

	if (g_pFileAgent_LR == NULL)
	{
		cout << "File output_agent_LR.csv cannot be opened." << endl;
		g_ProgramStop();
	}
	else
	{
		fprintf(g_pFileAgent_LR, "g_draw_seq_no, train_id, line_id, station_id, station_name, x_time_in_min, y_station_location\n");

		for (int a = 0; a < g_agent_vector.size(); a++)
		{
			CAgent* p_agent = &(g_agent_vector[a]);
			line_id = p_agent->agent_vector_seq_no + 1; //start from 1

			for (int fre = 0; fre < p_agent->frequency; fre++) //need to get frequency
			{
				g_draw_seq_no = 0;
				train_id += 1; //start from 1

				for (int l = 0; l < p_agent->g_output_link_NO_vector.size(); l++)
				{
					int link_no = p_agent->g_output_link_NO_vector[l];
					int from_node_no = g_link_vector[link_no].from_node_seq_no;
					int to_node_no = g_link_vector[link_no].to_node_seq_no;
					station_id = g_node_vector[from_node_no].station_id;
					station_name = g_node_vector[from_node_no].name;
					int temp_x_time_in_min = min(p_agent->g_output_link_TA_vector[l] + fre * ceil(g_cycle_length / p_agent->frequency), g_number_of_simulation_intervals - 1);
					int temp_x_to_node_time_in_min = min(p_agent->g_output_link_TD_vector[l] + fre * ceil(g_cycle_length / p_agent->frequency), g_number_of_simulation_intervals - 1);

					y_station_location = g_node_vector[from_node_no].x;
					temp_y_station_location = ceil(y_station_location);

					if (temp_x_time_in_min < g_cycle_length)
					{
						x_time_in_min = temp_x_time_in_min;
						fprintf(g_pFileAgent_LR, "%d,%d,%d,%d,%s,%d,%d\n", g_draw_seq_no, train_id, line_id, station_id, station_name.c_str(), x_time_in_min, temp_y_station_location);
					}
					else if (temp_x_time_in_min % g_cycle_length == 0 && temp_x_time_in_min > 0) //could be 1, 2,...
					{
						if (l == 0)
						{
							x_time_in_min = 0;
							fprintf(g_pFileAgent_LR, "%d,%d,%d,%d,%s,%d,%d\n", g_draw_seq_no, train_id, line_id, station_id, station_name.c_str(), x_time_in_min, temp_y_station_location);
						}
						else
						{
							x_time_in_min = g_cycle_length;
							fprintf(g_pFileAgent_LR, "%d,%d,%d,%d,%s,%d,%d\n", g_draw_seq_no, train_id, line_id, station_id, station_name.c_str(), x_time_in_min, temp_y_station_location);

							g_draw_seq_no += 1;

							x_time_in_min = 0;
							fprintf(g_pFileAgent_LR, "%d,%d,%d,%d,%s,%d,%d\n", g_draw_seq_no, train_id, line_id, station_id, station_name.c_str(), x_time_in_min, temp_y_station_location);
						}
					}
					else if (temp_x_time_in_min > g_cycle_length)
					{
						x_time_in_min = temp_x_time_in_min - (temp_x_time_in_min / g_cycle_length) * g_cycle_length;
						fprintf(g_pFileAgent_LR, "%d,%d,%d,%d,%s,%d,%d\n", g_draw_seq_no, train_id, line_id, station_id, station_name.c_str(), x_time_in_min, temp_y_station_location);
					}

					if ((temp_x_to_node_time_in_min / g_cycle_length) > (temp_x_time_in_min / g_cycle_length) && (temp_x_to_node_time_in_min % g_cycle_length) != 0) //cross the boundary 
					{
						//one node at the right boundary
						x_time_in_min = g_cycle_length;
						y_station_location = y_station_location + (g_node_vector[to_node_no].x - g_node_vector[from_node_no].x) * ((((temp_x_to_node_time_in_min / g_cycle_length) * g_cycle_length - temp_x_time_in_min) * 1.0f)
							/ ((temp_x_to_node_time_in_min - temp_x_time_in_min) * 1.0f));
						temp_y_station_location = ceil(y_station_location);
						fprintf(g_pFileAgent_LR, "%d,%d,%d,%d,%s,%d,%d\n", g_draw_seq_no, train_id, line_id, station_id, station_name.c_str(), x_time_in_min, temp_y_station_location);

						//another node at the left boundary
						g_draw_seq_no += 1;
						x_time_in_min = 0;
						fprintf(g_pFileAgent_LR, "%d,%d,%d,%d,%s,%d,%d\n", g_draw_seq_no, train_id, line_id, station_id, station_name.c_str(), x_time_in_min, temp_y_station_location);
					}

					if (l == (p_agent->g_output_link_NO_vector.size() - 1)) //add right node
					{
						station_id = g_node_vector[to_node_no].station_id;
						station_name = g_node_vector[to_node_no].name;

						if (temp_x_to_node_time_in_min > g_cycle_length)
						{
							x_time_in_min = temp_x_to_node_time_in_min - (temp_x_to_node_time_in_min / g_cycle_length) * g_cycle_length;

							if (x_time_in_min == 0)
							{
								x_time_in_min = g_cycle_length;
							}
						}
						else if (temp_x_to_node_time_in_min == g_cycle_length)
						{
							x_time_in_min = g_cycle_length;
						}
						else
						{
							x_time_in_min = temp_x_to_node_time_in_min;
						}

						y_station_location = g_node_vector[to_node_no].x;
						temp_y_station_location = ceil(y_station_location);
						fprintf(g_pFileAgent_LR, "%d,%d,%d,%d,%s,%d,%d\n", g_draw_seq_no, train_id, line_id, station_id, station_name.c_str(), x_time_in_min, temp_y_station_location);
					}
				}
			}
		}

		fclose(g_pFileAgent_LR);
	}
}

void g_OutAgentCSVFile_FromSimulation_ADMM()
{
	FILE* g_pFileAgent_ADMM = NULL;

	int train_id = 0;
	int line_id = 0;
	int station_id = 0;
	std::string station_name;
	int x_time_in_min = 0;
	float y_station_location = 0;
	int g_draw_seq_no = 0;
	int temp_y_station_location = 0;

	g_pFileAgent_ADMM = fopen("output_agent_ADMM.csv", "w");

	if (g_pFileAgent_ADMM == NULL)
	{
		cout << "File output_agent_ADMM.csv cannot be opened." << endl;
		g_ProgramStop();
	}
	else
	{
		fprintf(g_pFileAgent_ADMM, "g_draw_seq_no, train_id, line_id, station_id, station_name, x_time_in_min, y_station_location\n");

		for (int a = 0; a < g_agent_vector.size(); a++)
		{
			CAgent* p_agent = &(g_agent_vector[a]);
			line_id = p_agent->agent_vector_seq_no + 1; //start from 1

			for (int fre = 0; fre < p_agent->frequency; fre++) //need to get frequency
			{
				g_draw_seq_no = 0;
				train_id += 1; //start from 1

				for (int l = 0; l < p_agent->g_output_link_NO_vector.size(); l++)
				{
					int link_no = p_agent->g_output_link_NO_vector[l];
					int from_node_no = g_link_vector[link_no].from_node_seq_no;
					int to_node_no = g_link_vector[link_no].to_node_seq_no;
					station_id = g_node_vector[from_node_no].station_id;
					station_name = g_node_vector[from_node_no].name;
					int temp_x_time_in_min = min(p_agent->g_output_link_TA_vector[l] + fre * ceil(g_cycle_length / p_agent->frequency), g_number_of_simulation_intervals - 1);
					int temp_x_to_node_time_in_min = min(p_agent->g_output_link_TD_vector[l] + fre * ceil(g_cycle_length / p_agent->frequency), g_number_of_simulation_intervals - 1);

					y_station_location = g_node_vector[from_node_no].x;
					temp_y_station_location = ceil(y_station_location);

					if (temp_x_time_in_min < g_cycle_length)
					{
						x_time_in_min = temp_x_time_in_min;
						fprintf(g_pFileAgent_ADMM, "%d,%d,%d,%d,%s,%d,%d\n", g_draw_seq_no, train_id, line_id, station_id, station_name.c_str(), x_time_in_min, temp_y_station_location);
					}
					else if (temp_x_time_in_min % g_cycle_length == 0 && temp_x_time_in_min > 0) //could be 1, 2,...
					{
						if (l == 0)
						{
							x_time_in_min = 0;
							fprintf(g_pFileAgent_ADMM, "%d,%d,%d,%d,%s,%d,%d\n", g_draw_seq_no, train_id, line_id, station_id, station_name.c_str(), x_time_in_min, temp_y_station_location);
						}
						else
						{
							x_time_in_min = g_cycle_length;
							fprintf(g_pFileAgent_ADMM, "%d,%d,%d,%d,%s,%d,%d\n", g_draw_seq_no, train_id, line_id, station_id, station_name.c_str(), x_time_in_min, temp_y_station_location);

							g_draw_seq_no += 1;

							x_time_in_min = 0;
							fprintf(g_pFileAgent_ADMM, "%d,%d,%d,%d,%s,%d,%d\n", g_draw_seq_no, train_id, line_id, station_id, station_name.c_str(), x_time_in_min, temp_y_station_location);
						}
					}
					else if (temp_x_time_in_min > g_cycle_length)
					{
						x_time_in_min = temp_x_time_in_min - (temp_x_time_in_min / g_cycle_length) * g_cycle_length;
						fprintf(g_pFileAgent_ADMM, "%d,%d,%d,%d,%s,%d,%d\n", g_draw_seq_no, train_id, line_id, station_id, station_name.c_str(), x_time_in_min, temp_y_station_location);
					}

					if ((temp_x_to_node_time_in_min / g_cycle_length) > (temp_x_time_in_min / g_cycle_length) && (temp_x_to_node_time_in_min % g_cycle_length) != 0) //cross the boundary
					{
						//one node at the right boundary
						x_time_in_min = g_cycle_length;
						y_station_location = y_station_location + (g_node_vector[to_node_no].x - g_node_vector[from_node_no].x) * ((((temp_x_to_node_time_in_min / g_cycle_length) * g_cycle_length - temp_x_time_in_min) * 1.0f)
							/ ((temp_x_to_node_time_in_min - temp_x_time_in_min) * 1.0f));
						temp_y_station_location = ceil(y_station_location);
						fprintf(g_pFileAgent_ADMM, "%d,%d,%d,%d,%s,%d,%d\n", g_draw_seq_no, train_id, line_id, station_id, station_name.c_str(), x_time_in_min, temp_y_station_location);

						//another node at the left boundary
						g_draw_seq_no += 1;
						x_time_in_min = 0;
						fprintf(g_pFileAgent_ADMM, "%d,%d,%d,%d,%s,%d,%d\n", g_draw_seq_no, train_id, line_id, station_id, station_name.c_str(), x_time_in_min, temp_y_station_location);
					}

					if (l == (p_agent->g_output_link_NO_vector.size() - 1)) //add right node
					{
						station_id = g_node_vector[to_node_no].station_id;
						station_name = g_node_vector[to_node_no].name;

						if (temp_x_to_node_time_in_min > g_cycle_length)
						{
							x_time_in_min = temp_x_to_node_time_in_min - (temp_x_to_node_time_in_min / g_cycle_length) * g_cycle_length;

							if (x_time_in_min == 0)
							{
								x_time_in_min = g_cycle_length;
							}
						}
						else if (temp_x_to_node_time_in_min == g_cycle_length)
						{
							x_time_in_min = g_cycle_length;
						}
						else
						{
							x_time_in_min = temp_x_to_node_time_in_min;
						}

						y_station_location = g_node_vector[to_node_no].x;
						temp_y_station_location = ceil(y_station_location);
						fprintf(g_pFileAgent_ADMM, "%d,%d,%d,%d,%s,%d,%d\n", g_draw_seq_no, train_id, line_id, station_id, station_name.c_str(), x_time_in_min, temp_y_station_location);
					}
				}
			}
		}

		fclose(g_pFileAgent_ADMM);
	}
}

bool g_upper_bound_solution_feasibility_check(int g_upper_bound_solution_check_flag)  //1 check; 0 do not check
{
	if (g_upper_bound_solution_check_flag == 0)
	{
		return true;
	}

	//reset the link (i, j, t) visit counts and usage flag
	for (int l = 0; l < g_number_of_links; l++)
	{
		for (int t = 0; t < g_number_of_simulation_intervals; t++)
		{
			g_link_vector[l].g_link_time_visit_counts_departure_in_upper_bound[t] = 0;
			g_link_vector[l].g_link_time_visit_counts_arrival_in_upper_bound[t] = 0;

			for (int a = 0; a < g_number_of_agents; a++)
			{
				g_link_vector[l].g_link_time_train_visit_flag_departure_in_upper_bound[t][a] = 0;
				g_link_vector[l].g_link_time_train_visit_flag_arrival_in_upper_bound[t][a] = 0;
			}
		}
	}

	// update the link (i, j, t) visit counts 
	for (int a = 0; a < g_agent_vector.size(); a++)
	{
		CAgent* p_agent = &(g_agent_vector[a]);

		if (p_agent->path_new_link_id_vector_upper_bound.size() == 0)
		{
			return false;
		}

		for (int l = 0; l < p_agent->path_new_link_id_vector_upper_bound.size(); l++)  // for each link in the path of this agent
		{
			int link_seq_no = p_agent->path_new_link_id_vector_upper_bound[l];

			//mark the link travel time in the same direction
			int TA = p_agent->path_link_TA_vector_upper_bound[l];
			int TD = p_agent->path_link_TD_vector_upper_bound[l];
			int TA_left_time = TA;
			int TD_right_time = TD;
			int TA_left_time_headway = 0;
			int TD_right_time_headway = 0;
			int next_link_no = 0;
			int former_link_no = 0;

			if (g_link_vector[link_seq_no].link_type != siding_track_type && g_link_vector[link_seq_no].link_type != dummy_track_type) // 11 is the station track
			{
				if (l == 0) //first section
				{
					TA_left_time_headway = g_departure_headway_stop;
					next_link_no = p_agent->path_new_link_id_vector_upper_bound[l + 1];

					if (g_link_vector[next_link_no].link_type == siding_track_type) //a arrival event
					{
						TD_right_time_headway = g_arrival_headway_stop;
					}
					else
					{
						TD_right_time_headway = g_arrival_headway_passing;
					}
				}
				else if (l == p_agent->path_new_link_id_vector_upper_bound.size() - 1) //last section
				{
					TD_right_time_headway = g_arrival_headway_stop;

					former_link_no = p_agent->path_new_link_id_vector_upper_bound[l - 1];

					if (g_link_vector[former_link_no].link_type == siding_track_type) //a departure event
					{
						TA_left_time_headway = g_departure_headway_stop;
					}
					else //a passing event
					{
						TA_left_time_headway = g_departure_headway_passing;
					}
				}
				else //intermediate section
				{
					next_link_no = p_agent->path_new_link_id_vector_upper_bound[l + 1];

					if (g_link_vector[next_link_no].link_type == siding_track_type) //a arrival event
					{
						TD_right_time_headway = g_arrival_headway_stop;
					}
					else
					{
						TD_right_time_headway = g_arrival_headway_passing;
					}

					former_link_no = p_agent->path_new_link_id_vector_upper_bound[l - 1];

					if (g_link_vector[former_link_no].link_type == siding_track_type) //a departure event
					{
						TA_left_time_headway = g_departure_headway_stop;
					}
					else //a passing event
					{
						TA_left_time_headway = g_departure_headway_passing;
					}
				}

				for (int time = TA_left_time; time <= min(TA_left_time + TA_left_time_headway - 1, g_number_of_simulation_intervals - 1); time++)
				{
					if (g_link_vector[link_seq_no].g_link_time_train_visit_flag_departure_in_upper_bound[time][a] == 0)
					{
						for (int fre = 0; fre < p_agent->frequency; fre++)
						{
							for (int h = 0; h <= H; h++)
							{
								int temp_time = min(time + fre * ceil(g_cycle_length / p_agent->frequency) + h * g_cycle_length, g_number_of_simulation_intervals - 1);
								g_link_vector[link_seq_no].g_link_time_visit_counts_departure_in_upper_bound[temp_time] += 1;
								g_link_vector[link_seq_no].g_link_time_train_visit_flag_departure_in_upper_bound[temp_time][a] = 1;
							}
						}
					}
				}

				for (int time = TD_right_time; time <= min(TD_right_time + TD_right_time_headway - 1, g_number_of_simulation_intervals - 1); time++)
				{
					if (g_link_vector[link_seq_no].g_link_time_train_visit_flag_arrival_in_upper_bound[time][a] == 0)
					{
						for (int fre = 0; fre < p_agent->frequency; fre++)
						{
							for (int h = 0; h <= H; h++)
							{
								int temp_time = min(time + fre * ceil(g_cycle_length / p_agent->frequency) + h * g_cycle_length, g_number_of_simulation_intervals - 1);
								g_link_vector[link_seq_no].g_link_time_visit_counts_arrival_in_upper_bound[temp_time] += 1;
								g_link_vector[link_seq_no].g_link_time_train_visit_flag_arrival_in_upper_bound[temp_time][a] = 1;
							}
						}
					}
				}

			}
		}
	}

	bool bFeasibleFlag = true;

	for (int l = 0; l < g_link_vector.size(); l++)
	{
		if (g_link_vector[l].direction == 0)
		{
			if (g_link_vector[l].link_type != siding_track_type && g_link_vector[l].link_type != dummy_track_type)
			{
				for (int t = 0; t < g_number_of_simulation_intervals; t++)
				{
					int link_visit_counts = g_link_vector[l].g_link_time_visit_counts_departure_in_upper_bound[t];

					if (g_link_vector[l].same_link_id != 1000)
					{
						int same_link_no = g_internal_link_no_map[g_link_vector[l].same_link_id];
						link_visit_counts += g_link_vector[same_link_no].g_link_time_visit_counts_departure_in_upper_bound[t];
					}

					if (link_visit_counts > g_link_vector[l].time_depedent_capacity_matrix[t])
					{
						bFeasibleFlag = false;
						break;
					}

					link_visit_counts = g_link_vector[l].g_link_time_visit_counts_arrival_in_upper_bound[t];

					if (g_link_vector[l].same_link_id != 1000)
					{
						int same_link_no = g_internal_link_no_map[g_link_vector[l].same_link_id];
						link_visit_counts += g_link_vector[same_link_no].g_link_time_visit_counts_arrival_in_upper_bound[t];
					}

					if (link_visit_counts > g_link_vector[l].time_depedent_capacity_matrix[t])
					{
						bFeasibleFlag = false;
						break;
					}
				}

				if (bFeasibleFlag == false)
				{
					break;
				}

			}
		}
	}

	return bFeasibleFlag;
}

bool g_UpdateResourceUsageStatus()
{
	int w = 0;

	//reset the link (i, j, t) visit counts and usage flag
	for (int l = 0; l < g_number_of_links; l++)
	{
		for (int t = 0; t < g_number_of_simulation_intervals; t++)
		{
			g_link_vector[l].g_link_time_departure_visit_counts[t] = 0;
			g_link_vector[l].g_link_time_arrival_visit_counts[t] = 0;

			for (int a = 0; a < g_number_of_agents; a++)
			{
				g_link_vector[l].g_link_time_departure_train_visit_flag[t][a] = 0;
				g_link_vector[l].g_link_time_arrival_train_visit_flag[t][a] = 0;
			}
		}
	}

	// update the link (i, j, t) visit counts 
	for (int a = 0; a < g_agent_vector.size(); a++)
	{
		CAgent* p_agent = &(g_agent_vector[a]);

		for (int l = 0; l < p_agent->path_new_link_id_vector.size(); l++)  // for each link in the path of this agent
		{
			int TA_left_time = 0;
			int TD_right_time = 0;
			int TA_left_time_headway = 0;
			int TD_right_time_headway = 0;
			int next_link_no = 0;
			int former_link_no = 0;

			int link_seq_no = p_agent->path_new_link_id_vector[l];

			//mark the link travel time in the same direction
			int TA = p_agent->path_link_TA_vector[l];
			int TD = p_agent->path_link_TD_vector[l];
			TA_left_time = TA;
			TD_right_time = TD;

			if (g_link_vector[link_seq_no].link_type != siding_track_type && g_link_vector[link_seq_no].link_type != dummy_track_type)
			{
				if (l == 0) //first section
				{
					TA_left_time_headway = g_departure_headway_stop;
					next_link_no = p_agent->path_new_link_id_vector[l + 1];

					if (g_link_vector[next_link_no].link_type == siding_track_type) //a arrival event
					{
						TD_right_time_headway = g_arrival_headway_stop;
					}
					else
					{
						TD_right_time_headway = g_arrival_headway_passing;
					}
				}
				else if (l == p_agent->path_new_link_id_vector.size() - 1) //last section
				{
					TD_right_time_headway = g_arrival_headway_stop;

					former_link_no = p_agent->path_new_link_id_vector[l - 1];

					if (g_link_vector[former_link_no].link_type == siding_track_type) //a departure event
					{
						TA_left_time_headway = g_departure_headway_stop;
					}
					else //a passing event
					{
						TA_left_time_headway = g_departure_headway_passing;
					}
				}
				else //intermediate section
				{
					next_link_no = p_agent->path_new_link_id_vector[l + 1];

					if (g_link_vector[next_link_no].link_type == siding_track_type) //a arrival event
					{
						TD_right_time_headway = g_arrival_headway_stop;
					}
					else
					{
						TD_right_time_headway = g_arrival_headway_passing;
					}

					former_link_no = p_agent->path_new_link_id_vector[l - 1];

					if (g_link_vector[former_link_no].link_type == siding_track_type) //a departure event
					{
						TA_left_time_headway = g_departure_headway_stop;
					}
					else //a passing event
					{
						TA_left_time_headway = g_departure_headway_passing;
					}
				}

				for (int t = TA_left_time; t <= min(TA_left_time + TA_left_time_headway - 1, g_number_of_simulation_intervals - 1); t++)
				{
					if (g_link_vector[link_seq_no].g_link_time_departure_train_visit_flag[t][a] == 0)
					{
						for (int f = 0; f < p_agent->frequency; f++)
						{
							for (int h = 0; h <= H; h++)
							{
								int time = min(t + f * ceil(g_cycle_length / p_agent->frequency) + h * g_cycle_length, g_number_of_simulation_intervals - 1);
								g_link_vector[link_seq_no].g_link_time_departure_visit_counts[time] += 1;
								g_link_vector[link_seq_no].g_link_time_departure_train_visit_flag[time][a] = 1;
							}
						}
					}
				}

				for (int t = TD_right_time; t <= min(TD_right_time + TD_right_time_headway - 1, g_number_of_simulation_intervals - 1); t++)
				{
					if (g_link_vector[link_seq_no].g_link_time_arrival_train_visit_flag[t][a] == 0)
					{
						for (int f = 0; f < p_agent->frequency; f++)
						{
							for (int h = 0; h <= H; h++)
							{
								int time = min(t + f * ceil(g_cycle_length / p_agent->frequency) + h * g_cycle_length, g_number_of_simulation_intervals - 1);
								g_link_vector[link_seq_no].g_link_time_arrival_visit_counts[time] += 1;
								g_link_vector[link_seq_no].g_link_time_arrival_train_visit_flag[time][a] = 1;
							}
						}
					}
				}

			}
		}
	}

	return true;
}

void Output_SpaceTime_Path()
{
	if (g_LR_algorithm_debug_flag == 2)
	{
		for (int a = 0; a < g_agent_vector.size(); a++)
		{
			CAgent* p_agent = &(g_agent_vector[a]);

			for (int l = 0; l < p_agent->path_new_link_id_vector.size(); l++)  // for each link in the path of this agent
			{
				int link_seq_no = p_agent->path_new_link_id_vector[l];
				int from_node_no = g_link_vector[link_seq_no].from_node_seq_no;
				int to_node_no = g_link_vector[link_seq_no].to_node_seq_no;
				int from_node_id = g_internal_node_seq_no_to_node_id_map[from_node_no];
				int to_node_id = g_internal_node_seq_no_to_node_id_map[to_node_no];

				//mark the link travel time in the same direction
				int TA = p_agent->path_link_TA_vector[l];
				int TD = p_agent->path_link_TD_vector[l];

				fprintf(g_pFileDebugLog_LR, "agent_id = %d, link_seq_no = %d, from_node_id = %d, to_node_id = %d, TA = %d, TD = %d\n",
					p_agent->agent_id, link_seq_no, from_node_id, to_node_id, TA, TD);
			}
		}
	}
	else
	{
		return;
	}
}

void Output_SpaceTime_Path_Upper_Bound()
{
	if (g_LR_algorithm_debug_flag == 2)
	{
		for (int a = 0; a < g_agent_vector.size(); a++)
		{
			CAgent* p_agent = &(g_agent_vector[a]);

			for (int l = 0; l < p_agent->path_new_link_id_vector_upper_bound.size(); l++)  // for each link in the path of this agent
			{
				int link_seq_no = p_agent->path_new_link_id_vector_upper_bound[l];
				int from_node_no = g_link_vector[link_seq_no].from_node_seq_no;
				int to_node_no = g_link_vector[link_seq_no].to_node_seq_no;
				int from_node_id = g_internal_node_seq_no_to_node_id_map[from_node_no];
				int to_node_id = g_internal_node_seq_no_to_node_id_map[to_node_no];

				//mark the link travel time in the same direction
				int TA = p_agent->path_link_TA_vector_upper_bound[l];
				int TD = p_agent->path_link_TD_vector_upper_bound[l];

				fprintf(g_pFileDebugLog_LR, "agent_id = %d, link_seq_no = %d, from_node_id = %d, to_node_id = %d, TA = %d, TD = %d\n",
					p_agent->agent_id, link_seq_no, from_node_id, to_node_id, TA, TD);
			}
		}
	}
	else
	{
		return;
	}
}

void ADMM_Output_SpaceTime_Path_Upper_Bound()
{
	if (g_ADMM_algorithm_debug_flag == 2)
	{
		fprintf(g_pFileDebugLog_ADMM, "Space-time path of the trains in the upper bound solution:\n");

		for (int a = 0; a < g_agent_vector.size(); a++)
		{
			CAgent* p_agent = &(g_agent_vector[a]);

			for (int l = 0; l < p_agent->path_new_link_id_vector_upper_bound.size(); l++)  // for each link in the path of this agent
			{
				int link_seq_no = p_agent->path_new_link_id_vector_upper_bound[l];
				int from_node_no = g_link_vector[link_seq_no].from_node_seq_no;
				int to_node_no = g_link_vector[link_seq_no].to_node_seq_no;
				int from_node_id = g_internal_node_seq_no_to_node_id_map[from_node_no];
				int to_node_id = g_internal_node_seq_no_to_node_id_map[to_node_no];

				//mark the link travel time in the same direction
				int TA = p_agent->path_link_TA_vector_upper_bound[l];
				int TD = p_agent->path_link_TD_vector_upper_bound[l];

				fprintf(g_pFileDebugLog_ADMM, "agent_id = %d, link_seq_no = %d, from_node_id = %d, to_node_id = %d, TA = %d, TD = %d\n",
					p_agent->agent_id, link_seq_no, from_node_id, to_node_id, TA, TD);
			}
		}
	}
	else
	{
		return;
	}
}

void Output_SpaceTime_Path_ADMM()
{
	if (g_ADMM_algorithm_debug_flag == 2)
	{
		for (int a = 0; a < g_agent_vector.size(); a++)
		{
			CAgent* p_agent = &(g_agent_vector[a]);

			for (int l = 0; l < p_agent->path_new_link_id_vector.size(); l++)  // for each link in the path of this agent
			{
				int link_seq_no = p_agent->path_new_link_id_vector[l];
				int from_node_no = g_link_vector[link_seq_no].from_node_seq_no;
				int to_node_no = g_link_vector[link_seq_no].to_node_seq_no;
				int from_node_id = g_internal_node_seq_no_to_node_id_map[from_node_no];
				int to_node_id = g_internal_node_seq_no_to_node_id_map[to_node_no];

				//mark the link travel time in the same direction
				int TA = p_agent->path_link_TA_vector[l];
				int TD = p_agent->path_link_TD_vector[l];

				fprintf(g_pFileDebugLog_ADMM, "agent_id = %d, link_seq_no = %d, from_node_id = %d, to_node_id = %d, TA = %d, TD = %d\n",
					p_agent->agent_id, link_seq_no, from_node_id, to_node_id, TA, TD);
			}
		}
	}
	else
	{
		return;
	}
}

void g_UpperBound_ResourceMatrix_Initialization()
{
	for (int l = 0; l < g_link_vector.size(); l++)
	{
		g_link_vector[l].Setup_State_Dependent_Data_Matrix();
	}

	//initialization of state-time-space network
	CSTS_State element;

	g_STSStateVector.push_back(element);

	g_add_state_transition(0, 0, 0);

	int number_of_threads = 1;

	pSTSNetwork = new STSNetwork[number_of_threads];
	pSTSNetwork[0].AllocateSTSMemory(g_number_of_nodes, g_number_of_simulation_intervals, _MAX_STATES);
}

void g_UpperBound_ResourceMatrix_Clear()
{
	for (int l = 0; l < g_link_vector.size(); l++)
	{
		for (int t = 0; t < g_number_of_simulation_intervals; t++)
		{
			for (int s = 0; s < _MAX_STATES; s++)
			{
				if (g_link_vector[l].link_type == siding_track_type || g_link_vector[l].link_type == dummy_track_type)
				{
					g_link_vector[l].departure_state_dependent_travel_time_matrix[t][s] = g_link_vector[l].station_track_capacity - 1;
					g_link_vector[l].arrival_state_dependent_travel_time_matrix[t][s] = g_link_vector[l].station_track_capacity - 1;
				}
				else
				{
					g_link_vector[l].departure_state_dependent_travel_time_matrix[t][s] = 0;
					g_link_vector[l].arrival_state_dependent_travel_time_matrix[t][s] = 0;
				}
			}
		}
	}
}

void g_agent_vectors_clear()
{
	for (int a = 0; a < g_agent_vector.size(); a++)
	{
		g_agent_vector[a].path_new_link_id_vector.clear();
		g_agent_vector[a].path_link_TA_vector.clear();
		g_agent_vector[a].path_link_TD_vector.clear();
		g_agent_vector[a].path_new_node_id_vector.clear();
		g_agent_vector[a].path_new_timestamp_vector.clear();
		g_agent_vector[a].partial_schedule_node_vector.clear();
		g_agent_vector[a].partial_schedule_node_TD.clear();
		g_agent_vector[a].partial_schedule_node_TA.clear();
	}
}

bool lp_descending(int a, int b) //LP descending order
{
	if (g_agent_vector[a].m_DeviationRatio < g_agent_vector[b].m_DeviationRatio)
		return true;
	return false;
}

void g_train_ranking_by_LP_ADMM(int LR_iteration)
{
	g_agent_sequence_by_LP_ADMM.clear();

	for (int a = g_agent_vector.size() - 1; a >= 0; a--)
	{
		g_agent_sequence_by_LP_ADMM.push_back(a); //store the train ID first
	}

	for (int a = 0; a < g_agent_vector.size(); a++)
	{
		float free_trip_time = g_agent_vector[a].free_flow_travel_time;
		float actual_trip_time = g_agent_vector[a].travel_time_in_dual_solution;

		g_agent_vector[a].m_DeviationRatio = (actual_trip_time - free_trip_time) / free_trip_time;
	}

	sort(g_agent_sequence_by_LP_ADMM.begin(), g_agent_sequence_by_LP_ADMM.end(), lp_descending);
}

//deduce_flag: 0 LR; 1 priority rule; 2 ADMM
bool g_DeduceToProblemFeasibleSolution(STSNetwork pSTSNetwork[], std::vector<int> g_agent_sequences, int deduce_flag)
{
	int number_of_threads = 1;
	int thread_no = 0;

	//train sequences based on lagrangian profit

	for (int a = 0; a < g_agent_sequences.size(); a++)
	{
		int agent_no = g_agent_sequences[a];

		CAgent* p_agent = &(g_agent_vector[agent_no]);

		int internal_origin_node_seq_no = g_internal_node_seq_no_map[p_agent->origin_node_id];  // map external node number to internal node seq no. 

		pSTSNetwork[thread_no].m_agent_vector.clear();
		pSTSNetwork[thread_no].m_origin_node = internal_origin_node_seq_no;
		pSTSNetwork[thread_no].m_departure_time_beginning = p_agent->earliest_departure_time;

		pSTSNetwork[thread_no].m_arrival_time_ending = g_number_of_intervals_in_master_schedule - 1;
		pSTSNetwork[thread_no].m_agent_vector.push_back(agent_no);

		int LR_iteration = 0;

		pSTSNetwork[thread_no].find_STS_path_for_agents_assigned_for_this_thread(number_of_threads, LR_iteration, deduce_flag);
	}

	return true;
}

bool g_LR_Optimization(STSNetwork pSTSNetwork[])
{
	cout << "Lagrangian Relaxation Optimization..." << endl;

	g_SolutionStartTime = CTime::GetCurrentTime();

	LR_Initialization_start = clock();

	int number_of_threads = 1;
	int deduce_flag = 0;
	int w = 0;

	g_best_upper_bound = 99999;
	g_best_lower_bound = -99999;
	optimality_gap = 0;
	freStop_best_feas_Upper_bound_flag = 1;
	g_agent_vectors_clear();

	for (int l = 0; l < g_link_vector.size(); l++)
	{
		for (int t = 0; t < g_number_of_simulation_intervals; t++)
		{
			//initial value of LR multipliers
			g_link_vector[l].state_time_dependent_departure_LR_multiplier_matrix[t][w] = g_initial_LR_multiplier_value;
			g_link_vector[l].state_time_dependent_arrival_LR_multiplier_matrix[t][w] = g_initial_LR_multiplier_value;
		}
	}

	//reset the link (i, j, t) visit counts and usage flag
	for (int l = 0; l < g_number_of_links; l++)
	{
		for (int t = 0; t < g_number_of_simulation_intervals; t++)
		{
			g_link_vector[l].g_link_time_departure_visit_counts[t] = 0;
			g_link_vector[l].g_link_time_arrival_visit_counts[t] = 0;

			for (int a = 0; a < g_number_of_agents; a++)
			{
				g_link_vector[l].g_link_time_departure_train_visit_flag[t][a] = 0;
				g_link_vector[l].g_link_time_arrival_train_visit_flag[t][a] = 0;
			}
		}
	}

	if (g_LR_algorithm_debug_flag == 3)
	{
		fprintf(g_LR_iteration_Log, "iteration number, cpu time, best lower bound, best upper bound, optimality_gap\n");
		fprintf(g_LR_algorithmic_times_Log, "Initialization, Lower bound solution generation, Upper bound solution generation, Lagrangian multipliers updating\n");
	}

	LR_Initialization_end = clock();
	LR_Initialization_time = (LR_Initialization_end - LR_Initialization_start) / 1000.0f;

	//loop for each LR iteration
	for (int LR_iteration = 0; LR_iteration < g_number_of_LR_iterations; LR_iteration++)  // first loop
	{
		LR_LMU_start = clock();

		if ((LR_iteration + 1) % 200 == 0)
			cout << LR_iteration + 1 << "/" << g_number_of_LR_iterations << endl;

		g_stepSize = 1.0f / (LR_iteration + 1.0f);

		//keep the minimum step size
		if (g_stepSize < g_minimum_subgradient_step_size)
		{
			g_stepSize = g_minimum_subgradient_step_size;
		}

		for (int l = 0; l < g_link_vector.size(); l++)
		{
			if (g_link_vector[l].direction == 0)
			{
				for (int t = 0; t < g_number_of_simulation_intervals; t++)
				{
					//departure
					int link_visit_counts = g_link_vector[l].g_link_time_departure_visit_counts[t];

					if (g_link_vector[l].same_link_id != 1000)
					{
						int same_link_no = g_internal_link_no_map[g_link_vector[l].same_link_id];
						link_visit_counts += g_link_vector[same_link_no].g_link_time_departure_visit_counts[t];

						g_link_vector[l].state_time_dependent_departure_LR_multiplier_matrix[t][w] = max(0, g_link_vector[l].state_time_dependent_departure_LR_multiplier_matrix[t][w] +
							g_stepSize * (link_visit_counts - g_link_vector[l].time_depedent_capacity_matrix[t]));

						g_link_vector[same_link_no].state_time_dependent_departure_LR_multiplier_matrix[t][w] = g_link_vector[l].state_time_dependent_departure_LR_multiplier_matrix[t][w];
					}
					else
					{
						g_link_vector[l].state_time_dependent_departure_LR_multiplier_matrix[t][w] = max(0, g_link_vector[l].state_time_dependent_departure_LR_multiplier_matrix[t][w] +
							g_stepSize * (link_visit_counts - g_link_vector[l].time_depedent_capacity_matrix[t]));
					}

					if (g_LR_algorithm_debug_flag == 1)
					{
						fprintf(g_pFileDebugLog_LR, "link_no = %d, time = %d, link_visit_counts = %d, LR_multiplier = \t%0.2f\n", g_link_vector[l].link_seq_no, t,
							link_visit_counts, g_link_vector[l].state_time_dependent_departure_LR_multiplier_matrix[t][w]);
					}

					//arrival
					link_visit_counts = g_link_vector[l].g_link_time_arrival_visit_counts[t];

					if (g_link_vector[l].same_link_id != 1000)
					{
						int same_link_no = g_internal_link_no_map[g_link_vector[l].same_link_id];
						link_visit_counts += g_link_vector[same_link_no].g_link_time_arrival_visit_counts[t];

						g_link_vector[l].state_time_dependent_arrival_LR_multiplier_matrix[t][w] = max(0, g_link_vector[l].state_time_dependent_arrival_LR_multiplier_matrix[t][w] +
							g_stepSize * (link_visit_counts - g_link_vector[l].time_depedent_capacity_matrix[t]));

						g_link_vector[same_link_no].state_time_dependent_arrival_LR_multiplier_matrix[t][w] = g_link_vector[l].state_time_dependent_arrival_LR_multiplier_matrix[t][w];
					}
					else
					{
						g_link_vector[l].state_time_dependent_arrival_LR_multiplier_matrix[t][w] = max(0, g_link_vector[l].state_time_dependent_arrival_LR_multiplier_matrix[t][w] +
							g_stepSize * (link_visit_counts - g_link_vector[l].time_depedent_capacity_matrix[t]));
					}

					if (g_LR_algorithm_debug_flag == 1)
					{
						fprintf(g_pFileDebugLog_LR, "link_no = %d, time = %d, link_visit_counts = %d, LR_multiplier = \t%0.2f\n", g_link_vector[l].link_seq_no, t,
							link_visit_counts, g_link_vector[l].state_time_dependent_arrival_LR_multiplier_matrix[t][w]);
					}

				}
			}
		}

		LR_LMU_end = clock();
		LR_LMU_time += (LR_LMU_end - LR_LMU_start) / 1000.0f; //sum of the lr multipers updating time

		LR_LB_start = clock();

		// reset local LR lower bound
		float LR_global_lower_bound = 0;
		float total_price = 0;

		g_agent_vectors_clear();

		for (int a = 0; a < g_agent_vector.size(); a++)
		{
			int agent_no = g_agent_vector[a].agent_vector_seq_no;

			CAgent* p_agent = &(g_agent_vector[agent_no]);

			int internal_origin_node_seq_no = g_internal_node_seq_no_map[p_agent->origin_node_id];  // map external node number to internal node seq no. 
			int thread_no = 0;

			pSTSNetwork[thread_no].m_agent_vector.clear();
			pSTSNetwork[thread_no].m_origin_node = internal_origin_node_seq_no;
			pSTSNetwork[thread_no].m_departure_time_beginning = p_agent->earliest_departure_time;

			pSTSNetwork[thread_no].m_arrival_time_ending = g_number_of_intervals_in_master_schedule - 1;
			pSTSNetwork[thread_no].m_agent_vector.push_back(agent_no);

			float trip_price = pSTSNetwork[thread_no].find_STS_path_for_agents_assigned_for_this_thread_LR(number_of_threads, LR_iteration);

			total_price += trip_price;

			if (g_LR_algorithm_debug_flag == 2)
			{
				fprintf(g_pFileDebugLog_LR, "agent_id = %d, trip_price = %0.2f\n", agent_no, trip_price);
			}
		}

		float total_resource_price = 0;

		for (int l = 0; l < g_link_vector.size(); l++)
		{
			for (int t = 0; t < g_number_of_simulation_intervals; t++)
			{
				total_resource_price += g_link_vector[l].state_time_dependent_departure_LR_multiplier_matrix[t][w];
				total_resource_price += g_link_vector[l].state_time_dependent_arrival_LR_multiplier_matrix[t][w];
			}
		}

		LR_global_lower_bound = total_price - total_resource_price;
		g_best_lower_bound = max(g_best_lower_bound, LR_global_lower_bound);

		//check whether the feasible solution has been found by LR
		g_UpdateResourceUsageStatus();

		fprintf(g_pFileDebugLog_LR, "Space-time path of the trains in the lower bound solution:\n");
		Output_SpaceTime_Path();

		LR_LB_end = clock();
		LR_LB_time += (LR_LB_end - LR_LB_start) / 1000.0f; //sum of the lr lb generation time

		LR_UB_start = clock();

		bool is_success_in_finding_upperbound_solution = true;

		g_train_ranking_by_LP_ADMM(LR_iteration);

		g_agent_vectors_clear(); //clear the solutions in Lower bound
		g_UpperBound_ResourceMatrix_Clear(); //clear the upper bound resource matrix

		is_success_in_finding_upperbound_solution = g_DeduceToProblemFeasibleSolution(pSTSNetwork, g_agent_sequence_by_LP_ADMM, deduce_flag);

		bool is_upperbound_feasbile = true;
		int upper_bound_check_result = 1;
		is_upperbound_feasbile = g_upper_bound_solution_feasibility_check(g_upper_bound_solution_check_flag);

		float total_travel_time = 0;

		if (is_upperbound_feasbile == true)
		{
			for (int a = 0; a < g_agent_vector.size(); a++)
			{
				int n = g_agent_vector[a].path_new_link_id_vector_upper_bound.size() - 1;

				int actual_departure_time = g_agent_vector[a].path_link_TA_vector_upper_bound[0];
				int actual_arrival_time = g_agent_vector[a].path_link_TD_vector_upper_bound[n];

				total_travel_time += (actual_arrival_time - actual_departure_time) * g_agent_vector[a].frequency * 2;
			}

			g_best_upper_bound = min(g_best_upper_bound, total_travel_time);
		}

		int dummy_link_no = 0;
		int dummy_link_id = 0;

		vector<int>::iterator ret;

		for (int a = 0; a < g_agent_vector.size(); a++)
		{
			dummy_link_id = g_agent_vector[a].set_of_allowed_links_LR[g_agent_vector[a].set_of_allowed_links_LR.size() - 1];
			dummy_link_no = g_internal_link_no_map[dummy_link_id];

			ret = std::find(g_agent_vector[a].path_new_link_id_vector_upper_bound.begin(), g_agent_vector[a].path_new_link_id_vector_upper_bound.end(), dummy_link_no);
			if (ret != g_agent_vector[a].path_new_link_id_vector_upper_bound.end())
			{
				is_upperbound_feasbile = false;
				break;
			}
		}

		//only obtain the feasible Lagrangian relaxation paths
		if (is_upperbound_feasbile == true)
		{
			for (int a = 0; a < g_agent_vector.size(); a++)
			{
				dummy_link_id = g_agent_vector[a].set_of_allowed_links_LR[g_agent_vector[a].set_of_allowed_links_LR.size() - 1];
				dummy_link_no = g_internal_link_no_map[dummy_link_id];

				ret = std::find(g_agent_vector[a].path_new_link_id_vector_upper_bound.begin(), g_agent_vector[a].path_new_link_id_vector_upper_bound.end(), dummy_link_no);
				if (ret != g_agent_vector[a].path_new_link_id_vector_upper_bound.end())
				{
					freStop_best_feas_Upper_bound_flag = 0;
					break;
				}
			}

			if (total_travel_time <= g_best_upper_bound) //update the upper bound solution
			{
				for (int a = 0; a < g_agent_vector.size(); a++)
				{
					g_agent_vector[a].g_output_link_NO_vector = g_agent_vector[a].path_new_link_id_vector_upper_bound;
					g_agent_vector[a].g_output_link_TD_vector = g_agent_vector[a].path_link_TD_vector_upper_bound;
					g_agent_vector[a].g_output_link_TA_vector = g_agent_vector[a].path_link_TA_vector_upper_bound;
				}
			}
		}
		else
		{
			upper_bound_check_result = 0;
		}

		fprintf(g_pFileDebugLog_LR, "Space-time path of the trains in the upper bound solution:\n");
		Output_SpaceTime_Path_Upper_Bound();

		LR_UB_end = clock();
		LR_UB_time += (LR_UB_end - LR_UB_start) / 1000.0f; //sum of the lr ub generation time

		optimality_gap = (g_best_upper_bound - g_best_lower_bound) / g_best_lower_bound;

		if (g_output_log_flag == true)
		{
			fprintf(g_pFileDebugLog_LR, "LR_iteration = %d, stepSize = %0.2f, total_price = %0.2f, total_resource_price = %0.2f, g_best_lower_bound = %0.2f, g_best_upper_bound = %0.2f, optimality_gap = %0.4f, upper_bound_check_result = %d\n",
				LR_iteration, g_stepSize, total_price, total_resource_price, g_best_lower_bound / 2.0f, g_best_upper_bound / 2.0f, optimality_gap, upper_bound_check_result);
			fprintf(g_pFileDebugLog_LR, "------------------------------------------------------------------------------------------------------------------------\n");
		}

		if (g_LR_algorithm_debug_flag == 3)
		{
			CTimeSpan ctime = CTime::GetCurrentTime() - g_SolutionStartTime;
			fprintf(g_LR_iteration_Log, "%d, %lld, %0.2f, %0.2f, %0.4f\n", LR_iteration, ctime.GetTotalSeconds(), g_best_lower_bound / 2.0f, g_best_upper_bound / 2.0f, optimality_gap);
		}

		g_agent_vectors_clear(); //clear the solutions in Upper bound

	} //End for LR iteration

	if (g_LR_algorithm_debug_flag == 3)
	{
		fprintf(g_LR_algorithmic_times_Log, "%0.2f, %0.2f, %0.2f, %0.2f\n", LR_Initialization_time, LR_LB_time, LR_UB_time, LR_LMU_time);
	}

	g_OutAgentCSVFile_FromSimulation_LR();

	g_LR_FreStop_best_lower_bound.push_back(g_best_lower_bound);
	g_LR_FreStop_best_Upper_bound.push_back(g_best_upper_bound);
	g_LR_FreStop_best_Optimality_gap.push_back(optimality_gap);
	g_LR_FreStop_best_feas_Upper_bound_flag.push_back(freStop_best_feas_Upper_bound_flag);

	return true;
}

bool g_ADMM_Optimization(STSNetwork pSTSNetwork[])
{
	cout << "ADMM Optimization..." << endl;

	g_SolutionStartTime = CTime::GetCurrentTime();
	ADMM_Initialization_start = clock();

	int number_of_threads = 1;
	int w = 0;
	int deduce_flag = 2;

	g_best_upper_bound = 99999;
	g_best_lower_bound = -99999;
	optimality_gap = 0;
	freStop_best_feas_Lower_bound_flag = 1;
	g_best_ADMM_Feasible_lower_bound = 99999;
	int g_cumulative_unfeasible_iterations = 0;
	int g_cumulative_feasible_iterations = 0;

	g_agent_vectors_clear();

	for (int l = 0; l < g_link_vector.size(); l++)
	{
		for (int t = 0; t < g_number_of_simulation_intervals; t++)
		{
			for (int s = 0; s < _MAX_STATES; s++)
			{
				g_link_vector[l].state_time_dependent_departure_LR_multiplier_matrix[t][s] = g_initial_LR_multiplier_value;
				g_link_vector[l].state_time_dependent_arrival_LR_multiplier_matrix[t][s] = g_initial_LR_multiplier_value;

				g_link_vector[l].state_dependent_time_dependent_departure_ADMM_multiplier_matrix[t][s] = 0;
				g_link_vector[l].state_dependent_time_dependent_arrival_ADMM_multiplier_matrix[t][s] = 0;
			}
		}
	}

	//reset the link (i, j, t) visit counts and usage flag
	for (int l = 0; l < g_number_of_links; l++)
	{
		for (int t = 0; t < g_number_of_simulation_intervals; t++)
		{
			g_link_vector[l].g_link_time_departure_visit_counts[t] = 0;
			g_link_vector[l].g_link_time_arrival_visit_counts[t] = 0;

			for (int a = 0; a < g_number_of_agents; a++)
			{
				g_link_vector[l].g_link_time_departure_train_visit_flag[t][a] = 0;
				g_link_vector[l].g_link_time_arrival_train_visit_flag[t][a] = 0;
			}
		}
	}

	if (g_ADMM_algorithm_debug_flag == 3)
	{
		fprintf(g_ADMM_iteration_Log, "iteration number, cpu time, g_penalty_RHO_initial, g_primal_residual_error, penalty_term, best lower bound, best feasible lower bound, optimality_gap\n");

		fprintf(g_ADMM_algorithmic_times_Log, "Initialization, Lower bound solution generation, ADMM solution generation, Lagrangian multipliers and penalty parameter value updating\n");
	}

	//store the initial sequence
	vector<int> m_agent_no_vector;
	vector<int> temp_ADMM_optimization_sequence;
	for (int a = 0; a < g_agent_vector.size(); a++)
	{
		g_ADMM_optimization_sequence.push_back(a);
		m_agent_no_vector.push_back(a);
	}

	ADMM_Initialization_end = clock();
	ADMM_Initialization_time = (ADMM_Initialization_end - ADMM_Initialization_start) / 1000.0f;

	//loop for each ADMM iteration
	for (int ADMM_iteration = 0; ADMM_iteration < g_number_of_ADMM_iterations; ADMM_iteration++)  // first loop
	{
		ADMM_UB_start = clock();

		if ((ADMM_iteration + 1) % 200 == 0)
		{
			cout << ADMM_iteration + 1 << "/" << g_number_of_ADMM_iterations << endl;
		}

		// reset local LR lower bound
		float ADMM_total_travel_time_in_lower_bound = 0;
		float ADMM_global_lower_bound = 0;

		//g_agent_vectors_clear();
		unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();

		if (g_random_permutation == true)
		{
			if (ADMM_iteration == 0)
			{
				temp_ADMM_optimization_sequence = g_ADMM_optimization_sequence;
				g_ADMM_optimization_sequence = m_agent_no_vector;
			}
			else
			{
				temp_ADMM_optimization_sequence = g_ADMM_optimization_sequence;
				g_ADMM_optimization_sequence = m_agent_no_vector;
				std::shuffle(g_ADMM_optimization_sequence.begin(), g_ADMM_optimization_sequence.end(), std::default_random_engine(seed));
			}
		}

		for (int a = 0; a < g_ADMM_optimization_sequence.size(); a++)
		{
			int agent_no = g_ADMM_optimization_sequence[a];

			CAgent* p_agent = &(g_agent_vector[agent_no]);

			int internal_origin_node_seq_no = g_internal_node_seq_no_map[p_agent->origin_node_id];  // map external node number to internal node seq no. 
			int thread_no = 0;

			pSTSNetwork[thread_no].m_agent_vector.clear();
			pSTSNetwork[thread_no].m_origin_node = internal_origin_node_seq_no;
			pSTSNetwork[thread_no].m_departure_time_beginning = p_agent->earliest_departure_time;

			pSTSNetwork[thread_no].m_arrival_time_ending = g_number_of_intervals_in_master_schedule - 1;  // to do, use a predefined arrival time from input file
			pSTSNetwork[thread_no].m_agent_vector.push_back(agent_no);

			//update the ADMM multipiler according to k+1 before and k after
			int former_train_no = 0;
			vector<int> temp_agents_set;

			if (g_random_permutation == true)
			{
				if (a == 0) //first train in the sequence
				{
					former_train_no = temp_ADMM_optimization_sequence[g_ADMM_optimization_sequence.size() - 1]; //last train
				}
				else
				{
					former_train_no = g_ADMM_optimization_sequence[a - 1];
				}
			}
			else
			{
				if (a == 0) //first train in the sequence
				{
					former_train_no = g_ADMM_optimization_sequence[g_ADMM_optimization_sequence.size() - 1]; //last train
				}
				else
				{
					former_train_no = g_ADMM_optimization_sequence[a - 1];
				}
			}

			temp_agents_set.push_back(agent_no);
			temp_agents_set.push_back(former_train_no);

			for (int n = 0; n < temp_agents_set.size(); n++)
			{
				int agent = temp_agents_set[n];

				//reset the link (i, j, t) visit counts and usage flag
				for (int l = 0; l < g_number_of_links; l++)
				{
					for (int t = 0; t < g_number_of_simulation_intervals; t++)
					{
						g_link_vector[l].g_link_time_departure_train_visit_flag[t][agent] = 0;
						g_link_vector[l].g_link_time_arrival_train_visit_flag[t][agent] = 0;
					}
				}

				for (int l = 0; l < g_agent_vector[agent].path_new_link_id_vector.size(); l++)  // for each link in the path of this agent
				{
					int link_seq_no = g_agent_vector[agent].path_new_link_id_vector[l];

					//mark the link travel time in the same direction
					int TA = g_agent_vector[agent].path_link_TA_vector[l];
					int TD = g_agent_vector[agent].path_link_TD_vector[l];
					int TA_left_time = TA;
					int TD_right_time = TD;
					int TA_left_time_headway = 0;
					int TD_right_time_headway = 0;
					int next_link_no = 0;
					int former_link_no = 0;

					if (g_link_vector[link_seq_no].link_type != siding_track_type && g_link_vector[link_seq_no].link_type != dummy_track_type) // 11 is the station track
					{
						if (l == 0) //first section
						{
							TA_left_time_headway = g_departure_headway_stop;
							next_link_no = g_agent_vector[agent].path_new_link_id_vector[l + 1];

							if (g_link_vector[next_link_no].link_type == siding_track_type) //a arrival event
							{
								TD_right_time_headway = g_arrival_headway_stop;
							}
							else
							{
								TD_right_time_headway = g_arrival_headway_passing;
							}
						}
						else if (l == g_agent_vector[agent].path_new_link_id_vector.size() - 1) //last section
						{
							TD_right_time_headway = g_arrival_headway_stop;

							former_link_no = g_agent_vector[agent].path_new_link_id_vector[l - 1];

							if (g_link_vector[former_link_no].link_type == siding_track_type) //a departure event
							{
								TA_left_time_headway = g_departure_headway_stop;
							}
							else //a passing event
							{
								TA_left_time_headway = g_departure_headway_passing;
							}
						}
						else //intermediate section
						{
							next_link_no = g_agent_vector[agent].path_new_link_id_vector[l + 1];

							if (g_link_vector[next_link_no].link_type == siding_track_type) //a arrival event
							{
								TD_right_time_headway = g_arrival_headway_stop;
							}
							else
							{
								TD_right_time_headway = g_arrival_headway_passing;
							}

							former_link_no = g_agent_vector[agent].path_new_link_id_vector[l - 1];

							if (g_link_vector[former_link_no].link_type == siding_track_type) //a departure event
							{
								TA_left_time_headway = g_departure_headway_stop;
							}
							else //a passing event
							{
								TA_left_time_headway = g_departure_headway_passing;
							}
						}

						for (int t = TA_left_time; t <= min(TA_left_time + TA_left_time_headway - 1, g_number_of_simulation_intervals - 1); t++)
						{
							if (g_link_vector[link_seq_no].g_link_time_departure_train_visit_flag[t][agent] == 0)
							{
								for (int f = 0; f < g_agent_vector[agent].frequency; f++)
								{
									if (agent == agent_no) //minus
									{
										for (int h = 0; h <= H; h++)
										{
											int time = min(t + f * ceil(g_cycle_length / g_agent_vector[agent].frequency) + h * g_cycle_length, g_number_of_simulation_intervals - 1);
											g_link_vector[link_seq_no].g_link_time_departure_visit_counts[time] -= 1;
											g_link_vector[link_seq_no].g_link_time_departure_train_visit_flag[time][agent] = 1;
										}
									}
									else if (agent == former_train_no) //plus
									{
										for (int h = 0; h <= H; h++)
										{
											int time = min(t + f * ceil(g_cycle_length / g_agent_vector[agent].frequency) + h * g_cycle_length, g_number_of_simulation_intervals - 1);
											g_link_vector[link_seq_no].g_link_time_departure_visit_counts[time] += 1;
											g_link_vector[link_seq_no].g_link_time_departure_train_visit_flag[time][agent] = 1;
										}
									}
								}
							}
						}

						for (int t = TD_right_time; t <= min(TD_right_time + TD_right_time_headway - 1, g_number_of_simulation_intervals - 1); t++)
						{
							if (g_link_vector[link_seq_no].g_link_time_arrival_train_visit_flag[t][agent] == 0)
							{
								for (int f = 0; f < g_agent_vector[agent].frequency; f++)
								{
									if (agent == agent_no) //minus
									{
										for (int h = 0; h <= H; h++)
										{
											int time = min(t + f * ceil(g_cycle_length / g_agent_vector[agent].frequency) + h * g_cycle_length, g_number_of_simulation_intervals - 1);
											g_link_vector[link_seq_no].g_link_time_arrival_visit_counts[time] -= 1;
											g_link_vector[link_seq_no].g_link_time_arrival_train_visit_flag[time][agent] = 1;
										}
									}
									else if (agent == former_train_no) //plus
									{
										for (int h = 0; h <= H; h++)
										{
											int time = min(t + f * ceil(g_cycle_length / g_agent_vector[agent].frequency) + h * g_cycle_length, g_number_of_simulation_intervals - 1);
											g_link_vector[link_seq_no].g_link_time_arrival_visit_counts[time] += 1;
											g_link_vector[link_seq_no].g_link_time_arrival_train_visit_flag[time][agent] = 1;
										}
									}
								}
							}
						}

					}
				}
			}

			//update the values of ADMM multipliers
			for (int l = 0; l < g_link_vector.size(); l++)
			{
				if (g_link_vector[l].direction == 0)
				{
					for (int t = 0; t < g_number_of_simulation_intervals; t++)
					{
						//departure
						int link_visit_counts = g_link_vector[l].g_link_time_departure_visit_counts[t];

						if (g_link_vector[l].same_link_id != 1000)
						{
							int same_link_no = g_internal_link_no_map[g_link_vector[l].same_link_id];
							link_visit_counts += g_link_vector[same_link_no].g_link_time_departure_visit_counts[t];

							//g_link_vector[l].state_dependent_time_dependent_departure_ADMM_multiplier_matrix[t][w] = max(0, g_link_vector[l].state_time_dependent_departure_LR_multiplier_matrix[t][w] +
							//	(g_penalty_RHO_initial / 2.0f) * (2 * link_visit_counts - g_link_vector[l].time_depedent_capacity_matrix[t]));
							g_link_vector[l].state_dependent_time_dependent_departure_ADMM_multiplier_matrix[t][w] = g_link_vector[l].state_time_dependent_departure_LR_multiplier_matrix[t][w] +
								(g_penalty_RHO_initial / 2.0f) * max(0, (2 * link_visit_counts - g_link_vector[l].time_depedent_capacity_matrix[t]));

							g_link_vector[same_link_no].state_dependent_time_dependent_departure_ADMM_multiplier_matrix[t][w] = g_link_vector[l].state_dependent_time_dependent_departure_ADMM_multiplier_matrix[t][w];
						}
						else
						{
							//g_link_vector[l].state_dependent_time_dependent_departure_ADMM_multiplier_matrix[t][w] = max(0, g_link_vector[l].state_time_dependent_departure_LR_multiplier_matrix[t][w] +
								//(g_penalty_RHO_initial / 2.0f) * (2 * link_visit_counts - g_link_vector[l].time_depedent_capacity_matrix[t]));
							g_link_vector[l].state_dependent_time_dependent_departure_ADMM_multiplier_matrix[t][w] = g_link_vector[l].state_time_dependent_departure_LR_multiplier_matrix[t][w] +
								(g_penalty_RHO_initial / 2.0f) * max(0, (2 * link_visit_counts - g_link_vector[l].time_depedent_capacity_matrix[t]));

						}

						//arrival
						link_visit_counts = g_link_vector[l].g_link_time_arrival_visit_counts[t];

						if (g_link_vector[l].same_link_id != 1000)
						{
							int same_link_no = g_internal_link_no_map[g_link_vector[l].same_link_id];
							link_visit_counts += g_link_vector[same_link_no].g_link_time_arrival_visit_counts[t];

							//g_link_vector[l].state_dependent_time_dependent_arrival_ADMM_multiplier_matrix[t][w] = max(0, g_link_vector[l].state_time_dependent_arrival_LR_multiplier_matrix[t][w] +
							//	(g_penalty_RHO_initial / 2.0f) * (2 * link_visit_counts - g_link_vector[l].time_depedent_capacity_matrix[t]));
							g_link_vector[l].state_dependent_time_dependent_arrival_ADMM_multiplier_matrix[t][w] = g_link_vector[l].state_time_dependent_arrival_LR_multiplier_matrix[t][w] +
								(g_penalty_RHO_initial / 2.0f) * max(0, (2 * link_visit_counts - g_link_vector[l].time_depedent_capacity_matrix[t]));

							g_link_vector[same_link_no].state_dependent_time_dependent_arrival_ADMM_multiplier_matrix[t][w] = g_link_vector[l].state_dependent_time_dependent_arrival_ADMM_multiplier_matrix[t][w];
						}
						else
						{
							//g_link_vector[l].state_dependent_time_dependent_arrival_ADMM_multiplier_matrix[t][w] = max(0, g_link_vector[l].state_time_dependent_arrival_LR_multiplier_matrix[t][w] +
							//	(g_penalty_RHO_initial / 2.0f) * (2 * link_visit_counts - g_link_vector[l].time_depedent_capacity_matrix[t]));
							g_link_vector[l].state_dependent_time_dependent_arrival_ADMM_multiplier_matrix[t][w] = g_link_vector[l].state_time_dependent_arrival_LR_multiplier_matrix[t][w] +
								(g_penalty_RHO_initial / 2.0f) * max(0, (2 * link_visit_counts - g_link_vector[l].time_depedent_capacity_matrix[t]));
						}
					}
				}
			}

			float trip_price = pSTSNetwork[thread_no].find_STS_path_for_agents_assigned_for_this_thread_ADMM(number_of_threads, ADMM_iteration);

			if (g_ADMM_algorithm_debug_flag == 2)
			{
				fprintf(g_pFileDebugLog_ADMM, "agent_id = %d, trip_price = %0.2f\n", agent_no, trip_price);
			}
		}

		ADMM_UB_end = clock();
		ADMM_UB_time += (ADMM_UB_end - ADMM_UB_start) / 1000.0f; //sum of the lr ub generation time

		ADMM_LMU_start = clock();

		//get the link counts
		//reset the link (i, j, t) visit counts and usage flag
		for (int l = 0; l < g_number_of_links; l++)
		{
			for (int t = 0; t < g_number_of_simulation_intervals; t++)
			{
				g_link_vector[l].g_link_time_departure_visit_counts[t] = 0;
				g_link_vector[l].g_link_time_arrival_visit_counts[t] = 0;

				for (int a = 0; a < g_number_of_agents; a++)
				{
					g_link_vector[l].g_link_time_departure_train_visit_flag[t][a] = 0;
					g_link_vector[l].g_link_time_arrival_train_visit_flag[t][a] = 0;
				}
			}
		}

		//caculate the link time visit counts
		for (int a = 0; a < g_agent_vector.size(); a++)
		{
			int agent_no = g_agent_vector[a].agent_vector_seq_no;

			for (int l = 0; l < g_agent_vector[a].path_new_link_id_vector.size(); l++)  // for each link in the path of this agent
			{
				int link_seq_no = g_agent_vector[a].path_new_link_id_vector[l];

				//mark the link travel time in the same direction
				int TA = g_agent_vector[a].path_link_TA_vector[l];
				int TD = g_agent_vector[a].path_link_TD_vector[l];
				int TA_left_time = TA;
				int TD_right_time = TD;
				int TA_left_time_headway = 0;
				int TD_right_time_headway = 0;
				int next_link_no = 0;
				int former_link_no = 0;

				if (g_link_vector[link_seq_no].link_type != siding_track_type && g_link_vector[link_seq_no].link_type != dummy_track_type) // 11 is the station track
				{
					if (l == 0) //first section
					{
						TA_left_time_headway = g_departure_headway_stop;
						next_link_no = g_agent_vector[a].path_new_link_id_vector[l + 1];

						if (g_link_vector[next_link_no].link_type == siding_track_type) //a arrival event
						{
							TD_right_time_headway = g_arrival_headway_stop;
						}
						else
						{
							TD_right_time_headway = g_arrival_headway_passing;
						}
					}
					else if (l == g_agent_vector[a].path_new_link_id_vector.size() - 1) //last section
					{
						TD_right_time_headway = g_arrival_headway_stop;

						former_link_no = g_agent_vector[a].path_new_link_id_vector[l - 1];

						if (g_link_vector[former_link_no].link_type == siding_track_type) //a departure event
						{
							TA_left_time_headway = g_departure_headway_stop;
						}
						else //a passing event
						{
							TA_left_time_headway = g_departure_headway_passing;
						}
					}
					else //intermediate section
					{
						next_link_no = g_agent_vector[a].path_new_link_id_vector[l + 1];

						if (g_link_vector[next_link_no].link_type == siding_track_type) //a arrival event
						{
							TD_right_time_headway = g_arrival_headway_stop;
						}
						else
						{
							TD_right_time_headway = g_arrival_headway_passing;
						}

						former_link_no = g_agent_vector[a].path_new_link_id_vector[l - 1];

						if (g_link_vector[former_link_no].link_type == siding_track_type) //a departure event
						{
							TA_left_time_headway = g_departure_headway_stop;
						}
						else //a passing event
						{
							TA_left_time_headway = g_departure_headway_passing;
						}
					}

					for (int t = TA_left_time; t <= min(TA_left_time + TA_left_time_headway - 1, g_number_of_simulation_intervals - 1); t++)
					{
						if (g_link_vector[link_seq_no].g_link_time_departure_train_visit_flag[t][agent_no] == 0)
						{
							for (int f = 0; f < g_agent_vector[a].frequency; f++)
							{
								for (int h = 0; h <= H; h++)
								{
									int time = min(t + f * ceil(g_cycle_length / g_agent_vector[a].frequency) + h * g_cycle_length, g_number_of_simulation_intervals - 1);
									g_link_vector[link_seq_no].g_link_time_departure_visit_counts[time] += 1;
									g_link_vector[link_seq_no].g_link_time_departure_train_visit_flag[time][agent_no] = 1;
								}
							}
						}
					}

					for (int t = TD_right_time; t <= min(TD_right_time + TD_right_time_headway - 1, g_number_of_simulation_intervals - 1); t++)
					{
						if (g_link_vector[link_seq_no].g_link_time_arrival_train_visit_flag[t][agent_no] == 0)
						{
							for (int f = 0; f < g_agent_vector[a].frequency; f++)
							{
								for (int h = 0; h <= H; h++)
								{
									int time = min(t + f * ceil(g_cycle_length / g_agent_vector[a].frequency) + h * g_cycle_length, g_number_of_simulation_intervals - 1);
									g_link_vector[link_seq_no].g_link_time_arrival_visit_counts[time] += 1;
									g_link_vector[link_seq_no].g_link_time_arrival_train_visit_flag[time][agent_no] = 1;
								}
							}
						}
					}
				}
			}
		}

		//update the Lagrangian multiplers
		if (penalty_parameter_increasing_strategy == 0)
		{
			for (int l = 0; l < g_link_vector.size(); l++)
			{
				if (g_link_vector[l].direction == 0)
				{
					for (int t = 0; t < g_number_of_simulation_intervals; t++)
					{
						//departure
						int link_visit_counts = g_link_vector[l].g_link_time_departure_visit_counts[t];

						if (g_link_vector[l].same_link_id != 1000)
						{
							int same_link_no = g_internal_link_no_map[g_link_vector[l].same_link_id];
							link_visit_counts += g_link_vector[same_link_no].g_link_time_departure_visit_counts[t];

							g_link_vector[l].state_time_dependent_departure_LR_multiplier_matrix[t][w] = max(0, g_link_vector[l].state_time_dependent_departure_LR_multiplier_matrix[t][w] +
								g_gama_ratio * g_penalty_RHO_initial * (link_visit_counts - g_link_vector[l].time_depedent_capacity_matrix[t]));

							g_link_vector[same_link_no].state_time_dependent_departure_LR_multiplier_matrix[t][w] = g_link_vector[l].state_time_dependent_departure_LR_multiplier_matrix[t][w];
						}
						else
						{
							g_link_vector[l].state_time_dependent_departure_LR_multiplier_matrix[t][w] = max(0, g_link_vector[l].state_time_dependent_departure_LR_multiplier_matrix[t][w] +
								g_gama_ratio * g_penalty_RHO_initial * (link_visit_counts - g_link_vector[l].time_depedent_capacity_matrix[t]));
						}

						if (g_ADMM_algorithm_debug_flag == 1)
						{
							fprintf(g_pFileDebugLog_ADMM, "link_no = %d, time = %d, link_visit_counts = %d, LR_multiplier = \t%0.2f\n", g_link_vector[l].link_seq_no, t,
								link_visit_counts, g_link_vector[l].state_time_dependent_departure_LR_multiplier_matrix[t][w]);
						}

						//arrival
						link_visit_counts = g_link_vector[l].g_link_time_arrival_visit_counts[t];

						if (g_link_vector[l].same_link_id != 1000)
						{
							int same_link_no = g_internal_link_no_map[g_link_vector[l].same_link_id];
							link_visit_counts += g_link_vector[same_link_no].g_link_time_arrival_visit_counts[t];

							g_link_vector[l].state_time_dependent_arrival_LR_multiplier_matrix[t][w] = max(0, g_link_vector[l].state_time_dependent_arrival_LR_multiplier_matrix[t][w] +
								g_gama_ratio * g_penalty_RHO_initial * (link_visit_counts - g_link_vector[l].time_depedent_capacity_matrix[t]));

							g_link_vector[same_link_no].state_time_dependent_arrival_LR_multiplier_matrix[t][w] = g_link_vector[l].state_time_dependent_arrival_LR_multiplier_matrix[t][w];
						}
						else
						{
							g_link_vector[l].state_time_dependent_arrival_LR_multiplier_matrix[t][w] = max(0, g_link_vector[l].state_time_dependent_arrival_LR_multiplier_matrix[t][w] +
								g_gama_ratio * g_penalty_RHO_initial * (link_visit_counts - g_link_vector[l].time_depedent_capacity_matrix[t]));
						}

						if (g_ADMM_algorithm_debug_flag == 1)
						{
							fprintf(g_pFileDebugLog_ADMM, "link_no = %d, time = %d, link_visit_counts = %d, LR_multiplier = \t%0.2f\n", g_link_vector[l].link_seq_no, t,
								link_visit_counts, g_link_vector[l].state_time_dependent_arrival_LR_multiplier_matrix[t][w]);
						}
					}
				}
			}
		}
		else if (penalty_parameter_increasing_strategy == 1)
		{
			if (ADMM_iteration > number_of_lagrangian_iterations_for_ADMM)
			{
				for (int l = 0; l < g_link_vector.size(); l++)
				{
					if (g_link_vector[l].direction == 0)
					{
						for (int t = 0; t < g_number_of_simulation_intervals; t++)
						{
							//departure
							int link_visit_counts = g_link_vector[l].g_link_time_departure_visit_counts[t];

							if (g_link_vector[l].same_link_id != 1000)
							{
								int same_link_no = g_internal_link_no_map[g_link_vector[l].same_link_id];
								link_visit_counts += g_link_vector[same_link_no].g_link_time_departure_visit_counts[t];

								g_link_vector[l].state_time_dependent_departure_LR_multiplier_matrix[t][w] = max(0, g_link_vector[l].state_time_dependent_departure_LR_multiplier_matrix[t][w] +
									g_gama_ratio * g_penalty_RHO_initial * (link_visit_counts - g_link_vector[l].time_depedent_capacity_matrix[t]));

								g_link_vector[same_link_no].state_time_dependent_departure_LR_multiplier_matrix[t][w] = g_link_vector[l].state_time_dependent_departure_LR_multiplier_matrix[t][w];
							}
							else
							{
								g_link_vector[l].state_time_dependent_departure_LR_multiplier_matrix[t][w] = max(0, g_link_vector[l].state_time_dependent_departure_LR_multiplier_matrix[t][w] +
									g_gama_ratio * g_penalty_RHO_initial * (link_visit_counts - g_link_vector[l].time_depedent_capacity_matrix[t]));
							}

							if (g_ADMM_algorithm_debug_flag == 1)
							{
								fprintf(g_pFileDebugLog_ADMM, "link_no = %d, time = %d, link_visit_counts = %d, LR_multiplier = \t%0.2f\n", g_link_vector[l].link_seq_no, t,
									link_visit_counts, g_link_vector[l].state_time_dependent_departure_LR_multiplier_matrix[t][w]);
							}

							//arrival
							link_visit_counts = g_link_vector[l].g_link_time_arrival_visit_counts[t];

							if (g_link_vector[l].same_link_id != 1000)
							{
								int same_link_no = g_internal_link_no_map[g_link_vector[l].same_link_id];
								link_visit_counts += g_link_vector[same_link_no].g_link_time_arrival_visit_counts[t];

								g_link_vector[l].state_time_dependent_arrival_LR_multiplier_matrix[t][w] = max(0, g_link_vector[l].state_time_dependent_arrival_LR_multiplier_matrix[t][w] +
									g_gama_ratio * g_penalty_RHO_initial * (link_visit_counts - g_link_vector[l].time_depedent_capacity_matrix[t]));

								g_link_vector[same_link_no].state_time_dependent_arrival_LR_multiplier_matrix[t][w] = g_link_vector[l].state_time_dependent_arrival_LR_multiplier_matrix[t][w];
							}
							else
							{
								g_link_vector[l].state_time_dependent_arrival_LR_multiplier_matrix[t][w] = max(0, g_link_vector[l].state_time_dependent_arrival_LR_multiplier_matrix[t][w] +
									g_gama_ratio * g_penalty_RHO_initial * (link_visit_counts - g_link_vector[l].time_depedent_capacity_matrix[t]));
							}

							if (g_ADMM_algorithm_debug_flag == 1)
							{
								fprintf(g_pFileDebugLog_ADMM, "link_no = %d, time = %d, link_visit_counts = %d, LR_multiplier = \t%0.2f\n", g_link_vector[l].link_seq_no, t,
									link_visit_counts, g_link_vector[l].state_time_dependent_arrival_LR_multiplier_matrix[t][w]);
							}
						}
					}
				}
			}
			else
			{
				for (int l = 0; l < g_link_vector.size(); l++)
				{
					if (g_link_vector[l].direction == 0)
					{
						for (int t = 0; t < g_number_of_simulation_intervals; t++)
						{
							//departure
							int link_visit_counts = g_link_vector[l].g_link_time_departure_visit_counts[t];

							if (g_link_vector[l].same_link_id != 1000)
							{
								int same_link_no = g_internal_link_no_map[g_link_vector[l].same_link_id];
								link_visit_counts += g_link_vector[same_link_no].g_link_time_departure_visit_counts[t];

								g_link_vector[l].state_time_dependent_departure_LR_multiplier_matrix[t][w] = max(0, g_link_vector[l].state_time_dependent_departure_LR_multiplier_matrix[t][w] +
									g_stepSize * (link_visit_counts - g_link_vector[l].time_depedent_capacity_matrix[t]));

								g_link_vector[same_link_no].state_time_dependent_departure_LR_multiplier_matrix[t][w] = g_link_vector[l].state_time_dependent_departure_LR_multiplier_matrix[t][w];
							}
							else
							{
								g_link_vector[l].state_time_dependent_departure_LR_multiplier_matrix[t][w] = max(0, g_link_vector[l].state_time_dependent_departure_LR_multiplier_matrix[t][w] +
									g_stepSize * (link_visit_counts - g_link_vector[l].time_depedent_capacity_matrix[t]));
							}

							if (g_ADMM_algorithm_debug_flag == 1)
							{
								fprintf(g_pFileDebugLog_ADMM, "link_no = %d, time = %d, link_visit_counts = %d, LR_multiplier = \t%0.2f\n", g_link_vector[l].link_seq_no, t,
									link_visit_counts, g_link_vector[l].state_time_dependent_departure_LR_multiplier_matrix[t][w]);
							}

							//arrival
							link_visit_counts = g_link_vector[l].g_link_time_arrival_visit_counts[t];

							if (g_link_vector[l].same_link_id != 1000)
							{
								int same_link_no = g_internal_link_no_map[g_link_vector[l].same_link_id];
								link_visit_counts += g_link_vector[same_link_no].g_link_time_arrival_visit_counts[t];

								g_link_vector[l].state_time_dependent_arrival_LR_multiplier_matrix[t][w] = max(0, g_link_vector[l].state_time_dependent_arrival_LR_multiplier_matrix[t][w] +
									g_stepSize * (link_visit_counts - g_link_vector[l].time_depedent_capacity_matrix[t]));

								g_link_vector[same_link_no].state_time_dependent_arrival_LR_multiplier_matrix[t][w] = g_link_vector[l].state_time_dependent_arrival_LR_multiplier_matrix[t][w];
							}
							else
							{
								g_link_vector[l].state_time_dependent_arrival_LR_multiplier_matrix[t][w] = max(0, g_link_vector[l].state_time_dependent_arrival_LR_multiplier_matrix[t][w] +
									g_stepSize * (link_visit_counts - g_link_vector[l].time_depedent_capacity_matrix[t]));
							}

							if (g_ADMM_algorithm_debug_flag == 1)
							{
								fprintf(g_pFileDebugLog_ADMM, "link_no = %d, time = %d, link_visit_counts = %d, LR_multiplier = \t%0.2f\n", g_link_vector[l].link_seq_no, t,
									link_visit_counts, g_link_vector[l].state_time_dependent_arrival_LR_multiplier_matrix[t][w]);
							}
						}
					}
				}
			}
		}

		ADMM_LMU_end = clock();
		ADMM_LMU_time += (ADMM_LMU_end - ADMM_LMU_start) / 1000.0f; //sum of the lr multipers updating time

		ADMM_UB_start = clock();

		g_reoptimization_flag = false;
		int penalty_term = 0;

		//check whether the feasible solution has been found by ADMM
		bool bFeasibleFlag = true;

		for (int l = 0; l < g_link_vector.size(); l++)
		{
			if (g_link_vector[l].direction == 0)
			{
				for (int t = 0; t < g_number_of_simulation_intervals; t++)
				{
					//departure
					int link_visit_counts = g_link_vector[l].g_link_time_departure_visit_counts[t];

					if (g_link_vector[l].same_link_id != 1000)
					{
						int same_link_no = g_internal_link_no_map[g_link_vector[l].same_link_id];
						link_visit_counts += g_link_vector[same_link_no].g_link_time_departure_visit_counts[t];
					}

					if (link_visit_counts > g_link_vector[l].time_depedent_capacity_matrix[t])
					{
						bFeasibleFlag = false;
						penalty_term += 1;
					}

					//arrival
					link_visit_counts = g_link_vector[l].g_link_time_arrival_visit_counts[t];

					if (g_link_vector[l].same_link_id != 1000)
					{
						int same_link_no = g_internal_link_no_map[g_link_vector[l].same_link_id];
						link_visit_counts += g_link_vector[same_link_no].g_link_time_arrival_visit_counts[t];
					}

					if (link_visit_counts > g_link_vector[l].time_depedent_capacity_matrix[t])
					{
						bFeasibleFlag = false;
						penalty_term += 1;
					}
				}
			}
		}

		ADMM_UB_end = clock();
		ADMM_UB_time += (ADMM_UB_end - ADMM_UB_start) / 1000.0f; //sum of the lr ub generation time

		g_penalty_term_vector.push_back(penalty_term);

		g_primal_residual_error = min(g_primal_residual_error, penalty_term);

		//temporally store the path and time vectors
		for (int a = 0; a < g_agent_vector.size(); a++)
		{
			g_agent_vector[a].path_new_link_id_vector_temp.clear();
			g_agent_vector[a].path_link_TA_vector_temp.clear();
			g_agent_vector[a].path_link_TD_vector_temp.clear();

			g_agent_vector[a].path_new_link_id_vector_temp = g_agent_vector[a].path_new_link_id_vector;
			g_agent_vector[a].path_link_TA_vector_temp = g_agent_vector[a].path_link_TA_vector;
			g_agent_vector[a].path_link_TD_vector_temp = g_agent_vector[a].path_link_TD_vector;
		}

		ADMM_LB_start = clock();

		g_agent_vectors_clear();
		float total_price = 0;
		float total_resource_price = 0;

		for (int a = 0; a < g_agent_vector.size(); a++)
		{
			int agent_no = g_agent_vector[a].agent_vector_seq_no;

			CAgent* p_agent = &(g_agent_vector[agent_no]);

			int internal_origin_node_seq_no = g_internal_node_seq_no_map[p_agent->origin_node_id];  // map external node number to internal node seq no. 
			int thread_no = 0;

			pSTSNetwork[thread_no].m_agent_vector.clear();
			pSTSNetwork[thread_no].m_origin_node = internal_origin_node_seq_no;
			pSTSNetwork[thread_no].m_departure_time_beginning = p_agent->earliest_departure_time;

			pSTSNetwork[thread_no].m_arrival_time_ending = g_number_of_intervals_in_master_schedule - 1;  // to do, use a predefined arrival time from input file
			pSTSNetwork[thread_no].m_agent_vector.push_back(agent_no);

			float trip_price = pSTSNetwork[thread_no].find_STS_path_for_agents_assigned_for_this_thread_LR(number_of_threads, ADMM_iteration);

			total_price += trip_price;

			if (g_ADMM_algorithm_debug_flag == 2)
			{
				fprintf(g_pFileDebugLog_LR, "agent_id = %d, trip_price = %0.2f\n", agent_no, trip_price);
			}
		}

		for (int l = 0; l < g_link_vector.size(); l++)
		{
			for (int t = 0; t < g_number_of_simulation_intervals; t++)
			{
				total_resource_price += g_link_vector[l].state_time_dependent_departure_LR_multiplier_matrix[t][w];
				total_resource_price += g_link_vector[l].state_time_dependent_arrival_LR_multiplier_matrix[t][w];
			}
		}

		ADMM_global_lower_bound = total_price - total_resource_price;
		g_best_lower_bound = max(g_best_lower_bound, ADMM_global_lower_bound);

		ADMM_LB_end = clock();
		ADMM_LB_time += (ADMM_LB_end - ADMM_LB_start) / 1000.0f; //sum of the lr lb generation time

		//restore the path and time vectors
		for (int a = 0; a < g_agent_vector.size(); a++)
		{
			g_agent_vector[a].path_new_link_id_vector.clear();
			g_agent_vector[a].path_link_TA_vector.clear();
			g_agent_vector[a].path_link_TD_vector.clear();

			g_agent_vector[a].path_new_link_id_vector = g_agent_vector[a].path_new_link_id_vector_temp;
			g_agent_vector[a].path_link_TA_vector = g_agent_vector[a].path_link_TA_vector_temp;
			g_agent_vector[a].path_link_TD_vector = g_agent_vector[a].path_link_TD_vector_temp;
		}

		for (int a = 0; a < g_agent_vector.size(); a++)
		{
			int n = g_agent_vector[a].path_link_TD_vector.size() - 1;
			ADMM_total_travel_time_in_lower_bound += (g_agent_vector[a].path_link_TD_vector[n] - g_agent_vector[a].path_link_TA_vector[0])
				* g_agent_vector[a].frequency * 2;
		}

		if (bFeasibleFlag == true)
		{
			g_best_ADMM_Feasible_lower_bound = min(g_best_ADMM_Feasible_lower_bound, ADMM_total_travel_time_in_lower_bound);
		}

		optimality_gap = (g_best_ADMM_Feasible_lower_bound - g_best_lower_bound) / g_best_lower_bound;

		if (g_ADMM_algorithm_debug_flag == 3)
		{
			CTimeSpan ctime = CTime::GetCurrentTime() - g_SolutionStartTime;
			fprintf(g_ADMM_iteration_Log, "%d, %lld, %0.2f, %d, %d, %0.2f, %0.2f, %0.4f\n", ADMM_iteration, ctime.GetTotalSeconds(), g_penalty_RHO_initial, g_primal_residual_error, penalty_term, g_best_lower_bound / 2.0f, g_best_ADMM_Feasible_lower_bound / 2.0f, optimality_gap);
		}

		int dummy_link_no = 0;
		int dummy_link_id = 0;
		vector<int>::iterator ret;

		for (int a = 0; a < g_agent_vector.size(); a++)
		{
			dummy_link_id = g_agent_vector[a].set_of_allowed_links_LR[g_agent_vector[a].set_of_allowed_links_LR.size() - 1];
			dummy_link_no = g_internal_link_no_map[dummy_link_id];

			ret = std::find(g_agent_vector[a].path_link_seq_no_vector.begin(), g_agent_vector[a].path_link_seq_no_vector.end(), dummy_link_no);
			if (ret != g_agent_vector[a].path_link_seq_no_vector.end())
			{
				bFeasibleFlag = false;
				break;
			}
		}

		//only obtain the feasible ADMM paths
		if (bFeasibleFlag == true)
		{
			g_best_ADMM_Feasible_lower_bound = min(g_best_ADMM_Feasible_lower_bound, ADMM_total_travel_time_in_lower_bound);

			if (ADMM_total_travel_time_in_lower_bound <= g_best_ADMM_Feasible_lower_bound)
			{
				for (int a = 0; a < g_agent_vector.size(); a++)
				{
					g_agent_vector[a].g_output_link_NO_vector = g_agent_vector[a].path_new_link_id_vector;
					g_agent_vector[a].g_output_link_TD_vector = g_agent_vector[a].path_link_TD_vector;
					g_agent_vector[a].g_output_link_TA_vector = g_agent_vector[a].path_link_TA_vector;
				}
			}

			g_cumulative_feasible_iterations += 1;
		}
		else
		{
			g_cumulative_unfeasible_iterations += 1;
		}

		if (g_ADMM_feasible_and_stop_flag == 1) //need to regenerate sequences
		{
			if (penalty_parameter_increasing_strategy == 0)
			{
				if ((bFeasibleFlag == false) && (g_cumulative_unfeasible_iterations + 1) % g_penalty_RHO_update_iteration == 0)
				{
					g_reoptimization_flag = true;
					g_penalty_RHO_initial = g_penalty_RHO_initial + g_penalty_RHO_incremental;
				}

				if ((bFeasibleFlag == true) && (g_cumulative_feasible_iterations + 1) % g_penalty_RHO_update_iteration == 0)
				{
					g_reoptimization_flag = true;
					g_penalty_RHO_initial = g_penalty_RHO_initial + g_penalty_RHO_incremental;
				}
			}
			else if (penalty_parameter_increasing_strategy == 1)
			{
				if (ADMM_iteration > number_of_lagrangian_iterations_for_ADMM)
				{
					float temp_number1 = g_penalty_term_vector[ADMM_iteration] * g_penalty_term_vector[ADMM_iteration];
					float temp_number2 = g_penalty_gamma_ratio * (g_penalty_term_vector[ADMM_iteration - 1] * g_penalty_term_vector[ADMM_iteration - 1]);

					if (temp_number1 >= temp_number2)
					{
						g_reoptimization_flag = true;
						g_penalty_RHO_initial = g_penalty_RHO_initial + g_penalty_RHO_incremental;
					}
				}
				else
				{
					g_stepSize = 1.0f / (ADMM_iteration + 1.0f);

					//keep the minimum step size
					if (g_stepSize < g_minimum_subgradient_step_size)
					{
						g_stepSize = g_minimum_subgradient_step_size;
					}
				}
			}
		}

		fprintf(g_pFileDebugLog_ADMM, "Space-time path of the trains in the lower bound solution:\n");
		Output_SpaceTime_Path_ADMM();

		if (g_ADMM_feasible_and_stop_flag == 1) //need to regenerate sequences
		{
			if (bFeasibleFlag == true)
			{
				g_ADMM_optimization_sequence.clear();

				for (int a = g_agent_vector.size() - 1; a >= 0; a--)
				{
					g_ADMM_optimization_sequence.push_back(a); //store the train ID first
				}

				for (int a = 0; a < g_agent_vector.size(); a++)
				{
					float free_trip_time = g_agent_vector[a].free_flow_travel_time;
					float actual_trip_time = g_agent_vector[a].travel_time_in_dual_solution;

					g_agent_vector[a].m_DeviationRatio = (actual_trip_time - free_trip_time) / free_trip_time;
				}

				sort(g_ADMM_optimization_sequence.begin(), g_ADMM_optimization_sequence.end(), lp_descending);
			}
		}

		bool is_success_in_finding_upperbound_solution = true;
		bool is_upperbound_feasbile = false;
		int upper_bound_check_result = 0;

		if (g_deduce_feasible_upperbound_in_ADMM_flag == 1)
		{
			g_train_ranking_by_LP_ADMM(ADMM_iteration);

			//clear the upper bound resource matrix
			g_UpperBound_ResourceMatrix_Clear();

			is_success_in_finding_upperbound_solution = g_DeduceToProblemFeasibleSolution(pSTSNetwork, g_agent_sequence_by_LP_ADMM, deduce_flag);
			is_upperbound_feasbile = g_upper_bound_solution_feasibility_check(g_upper_bound_solution_check_flag);

			float total_travel_time = 0;

			for (int a = 0; a < g_agent_vector.size(); a++)
			{
				int n = g_agent_vector[a].path_link_TD_vector_upper_bound.size() - 1;

				int actual_departure_time = g_agent_vector[a].path_link_TA_vector_upper_bound[0];
				int actual_arrival_time = g_agent_vector[a].path_link_TD_vector_upper_bound[n];

				total_travel_time += (actual_arrival_time - actual_departure_time) * g_agent_vector[a].frequency * 2;
			}

			if (is_upperbound_feasbile == true)
			{
				g_best_upper_bound = min(g_best_upper_bound, total_travel_time);
				upper_bound_check_result = 1;
			}
			else
			{
				upper_bound_check_result = 0;
			}
		}

		int lower_bound_check_result = 1;

		if (bFeasibleFlag == true)
		{
			freStop_best_feas_Lower_bound_flag = 1;
		}
		else
		{
			lower_bound_check_result = 0;
		}

		ADMM_Output_SpaceTime_Path_Upper_Bound();

		if (g_output_log_flag == true)
		{
			fprintf(g_pFileDebugLog_ADMM, "ADMM_iteration = %d, g_best_lower_bound = %0.2f, primal_residual_error = %d, penalty_term = %d, g_penalty_RHO_initial = %0.2f, g_best_ADMM_Feasible_lower_bound = %0.2f, optimality_gap = %0.4f, lower_bound_check_result = %d\n",
				ADMM_iteration, g_best_lower_bound / 2.0f, g_primal_residual_error, penalty_term, g_penalty_RHO_initial, g_best_ADMM_Feasible_lower_bound / 2.0f, optimality_gap, lower_bound_check_result);

			fprintf(g_pFileDebugLog_ADMM, "------------------------------------------------------------------------------------------------------------------------\n");
		}

	} //End for ADMM iteration

	if (g_ADMM_algorithm_debug_flag == 3)
	{
		fprintf(g_ADMM_algorithmic_times_Log, "%0.2f, %0.2f, %0.2f, %0.2f\n", ADMM_Initialization_time, ADMM_LB_time, ADMM_UB_time, ADMM_LMU_time);
	}

	g_OutAgentCSVFile_FromSimulation_ADMM();

	g_ADMM_FreStop_best_lower_bound.push_back(g_best_lower_bound);
	g_ADMM_FreStop_best_Upper_bound.push_back(g_best_upper_bound);
	g_ADMM_FreStop_best_Optimality_gap.push_back(optimality_gap);
	g_ADMM_FreStop_best_feas_Lower_bound_flag.push_back(freStop_best_feas_Lower_bound_flag);

	return true;
}

bool g_periodical_timetable_optimization()
{
	clock_t start_LR, end_LR, total_LR;
	clock_t start_ADMM, end_ADMM, total_ADMM;

	//Lagrangian optimization for Lower and Upper bound
	start_LR = clock();

	g_LR_Optimization(pSTSNetwork);

	fprintf(g_pFileDebugLog_LR, "------------------------------------------------------------------------------------------------------------------------\n");
	cout << "End of LR Optimization Process. " << endl;

	end_LR = clock();
	total_LR = (end_LR - start_LR);

	//ADMM optimization for Lower and Upper bound
	start_ADMM = clock();

	g_ADMM_Optimization(pSTSNetwork);

	fprintf(g_pFileDebugLog_ADMM, "------------------------------------------------------------------------------------------------------------------------\n");

	cout << "End of ADMM Optimization Process. " << endl;

	end_ADMM = clock();
	total_ADMM = (end_ADMM - start_ADMM);

	cout << "CPU Running Time of LR = " << total_LR / 1000.0f << " seconds" << endl;
	cout << "CPU Running Time of ADMM = " << total_ADMM / 1000.0f << " seconds" << endl;

	fprintf(g_pFileOutputLog, "CPU Running Time of LR =,%.2f, seconds\n", total_LR / 1000.0f);
	fprintf(g_pFileOutputLog, "CPU Running Time of ADMM =,%.2f, seconds\n", total_ADMM / 1000.0f);

	return true;
}

int main(int argc, TCHAR* argv[], TCHAR* envp[])
{
	// definte timestamps
	clock_t start_t, end_t, total_t;
	clock_t start_Priority, end_Priority, total_Priority;
	int i;

	start_t = clock();
	g_pFileDebugLog = fopen("Debug.txt", "w");

	if (g_pFileDebugLog == NULL)
	{
		cout << "File Debug.txt cannot be opened." << endl;
		g_ProgramStop();
	}

	g_pFileDebugLog_LR = fopen("DebugLR.txt", "w");

	if (g_pFileDebugLog_LR == NULL)
	{
		cout << "File DebugLR.txt cannot be opened." << endl;
		g_ProgramStop();
	}

	g_pFileDebugLog_ADMM = fopen("DebugADMM.txt", "w");

	if (g_pFileDebugLog_ADMM == NULL)
	{
		cout << "File DebugADMM.txt cannot be opened." << endl;
		g_ProgramStop();
	}

	g_pFileOutputLog = fopen("output_solution.csv", "w");
	if (g_pFileOutputLog == NULL)
	{
		cout << "File output_solution.csv cannot be opened." << endl;
		g_ProgramStop();
	}

	g_LR_iteration_Log = fopen("output_LR_iteration_Log.csv", "w");

	if (g_LR_iteration_Log == NULL)
	{
		cout << "File output_LR_iteration_Log.csv cannot be opened." << endl;
		g_ProgramStop();
	}

	g_ADMM_iteration_Log = fopen("output_ADMM_iteration_Log.csv", "w");

	if (g_ADMM_iteration_Log == NULL)
	{
		cout << "File output_ADMM_iteration_Log.csv cannot be opened." << endl;
		g_ProgramStop();
	}

	g_LR_algorithmic_times_Log = fopen("output_LR_times_Log.csv", "w");

	if (g_LR_algorithmic_times_Log == NULL)
	{
		cout << "File output_LR_times_Log.csv cannot be opened." << endl;
		g_ProgramStop();
	}

	g_ADMM_algorithmic_times_Log = fopen("output_ADMM_times_Log.csv", "w");

	if (g_ADMM_algorithmic_times_Log == NULL)
	{
		cout << "File output_ADMM_times_Log.csv cannot be opened." << endl;
		g_ProgramStop();
	}

	// step 1: read input data of network and demand agent
	g_ReadInputData();

	// step 2: initialize the resource matrix in upper bound
	g_UpperBound_ResourceMatrix_Initialization();

	// step 3: Lagrangian and ADMM optimization for Lower and Upper bound
	g_periodical_timetable_optimization();

	end_t = clock();
	total_t = (end_t - start_t);

	cout << "Total CPU Running Time = " << total_t / 1000.0f << " seconds" << endl;
	fprintf(g_pFileOutputLog, "Total CPU Running Time =,%.2f, seconds\n", total_t / 1000.0);

	fclose(g_pFileOutputLog);
	fclose(g_pFileDebugLog);
	fclose(g_pFileDebugLog_LR);
	fclose(g_LR_iteration_Log);
	fclose(g_LR_algorithmic_times_Log);	
	fclose(g_ADMM_algorithmic_times_Log);
	
	cout << "End of Optimization " << endl;
	cout << "free memory.." << endl;
	cout << "done." << endl;

	g_node_vector.clear();
	g_link_vector.clear();

	getchar();

	return 1;
}
