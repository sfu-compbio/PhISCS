/*******************************************************************************
* Author: Ehsan Haghshenas
* Last update: Oct 19, 2017
*******************************************************************************/

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <set>
#include <getopt.h>
#include <unistd.h>
#include <sys/time.h>
#include <sys/resource.h>

using namespace std;

#define MAX_CELL 300
#define MAX_MUT 200

string  par_inputFile = "";
string  par_outDir = "";
double  par_fnRate = -1;
double  par_fnWeight;
double  par_fnWeight_neg;
double  par_fpRate = -1;
double  par_fpWeight;
double  par_fpWeight_neg;
double  par_const = 1000000;
double  par_precisionFactor = 1000;
int     par_colWeight = -1;
int     par_maxColRemove = 0;
string  par_bulkFile = "";
double  par_delta = 0.01;
string  par_maxSolver = "openwbo";
int     par_threads = 1;
bool    par_isTrueVAF = false;
bool    IS_PWCNF = true;
bool    INT_WEIGHTS = false;
string  MAXSAT_EXE;

int mat[MAX_CELL][MAX_MUT]; // the 0/1 matrix
vector<string> cellId;
vector<string> mutId;
int var_x[MAX_CELL][MAX_MUT]; // X variables for the maxSAT
int var_y[MAX_CELL][MAX_MUT]; // Y variables for the maxSAT; if(Iij==0) Yij=Xij and if(Iij==1) Yij=~Xij
int weight_x[MAX_CELL][MAX_MUT]; // weight of X variables
int var_b[MAX_MUT][MAX_MUT][2][2];
int var_k[MAX_MUT];
int var_a[MAX_MUT][MAX_MUT];
double vaf[MAX_MUT];
int vafP[MAX_MUT][MAX_MUT];
int vafT[MAX_MUT][MAX_MUT][MAX_MUT];
pair<int, int> map_y2ij[MAX_CELL * MAX_MUT + 10]; // maps Y variables to matrix position (row and column)
// pair<int, int> map_a2pq[10000]; // maps a variables to matrix position (row and column)
vector<string> clauseSoft; // the set of soft clauses for wcnf formulation
vector<string> clauseHard; // the set of soft clauses for wcnf formulation

int numMut = 0; // actual number of mutations (columns)
int numCell = 0; // actual number of cells (rows)
int numVarY = 0; // number of Y variables
int numVarX = 0; // number of X variables
int numVarB = 0; // number of B variables
int numVarK = 0; // number of K variables
int numVarA = 0; // number of A variables
int numVarW = 0; // number of W variables
int numZero = 0; // number of zeros in the input matrix
int numOne = 0; // number of ones in the input matrix
int numTwo = 0; // number of twos in the input matrix

// #define startVarY (0)
#define startVarX (numVarY)
#define startVarB (numVarY + numVarX)
#define startVarK (numVarY + numVarX + numVarB)
#define startVarA (numVarY + numVarX + numVarB + numVarK)
#define startVarW (numVarY + numVarX + numVarB + numVarK + numVarA)

string int2str(int n)
{
    ostringstream sout;
    sout<< n;
    return sout.str();
}
string double2str(double n)
{
    ostringstream sout;
    sout<< n;
    return sout.str();
}

int str2int(string s)
{
    int retVal;
    istringstream sin(s.c_str());
    sin >> retVal;
    return retVal;
}

double str2double(string s)
{
    double retVal;
    istringstream sin(s.c_str());
    sin >> retVal;
    return retVal;
}

inline double log10(double n)
{
    return log(n)/log(10);
}

// double getCpuTime()
// {
//  struct rusage t;
//  getrusage(RUSAGE_SELF, &t);
//  return t.ru_utime.tv_sec + t.ru_utime.tv_usec / 1000000.0 + t.ru_stime.tv_sec + t.ru_stime.tv_usec / 1000000.0;
// }

double getRealTime()
{
    struct timeval t;
    struct timezone tz;
    gettimeofday(&t, &tz);
    return t.tv_sec + t.tv_usec / 1000000.0;
}

string get_file_name(string path, bool removExtension = false)
{
    string fileName;
    size_t pos;
    // extract file name
    pos = path.find_last_of("/");
    if(pos != string::npos)
        fileName = path.substr(pos+1);
    else
        fileName = path;
    // remove extension
    if(removExtension)
    {
        pos = fileName.find_last_of(".");
        if(pos != string::npos)
            fileName = fileName.substr(0, pos);
    }
    return fileName;
}

string get_dir_path(string path)
{
    size_t pos = path.find_last_of("/");
    if(pos != string::npos)
    {
        return path.substr(0, pos);
    }
    else
    {
        return "";
    }
}

string get_exe_path()
{
  char path[10000];
  ssize_t count = readlink( "/proc/self/exe", path, 10000);
  return string(path, (count > 0) ? count : 0);
}

void print_usage()
{
    cout<< endl
        << "usage: csp_maxsat [-h] -f FILE -n FNRATE -p FPRATE -o OUTDIR" << endl
        << "                  [-c CONST] [-s SOLVER]" << endl;
        // << "                  [-m MAXMUT] [-b BULK] [-e DELTA] [-v] [-t THREADS]" << endl;
}

void print_help()
{
    cout<< endl
        << "Required arguments:" << endl
        << "   -f, --file     STR        Input matrix file" << endl
        << "   -n, --fnRate   FLT        The estimate for false negative rate" << endl
        << "   -p, --fpRate   FLT        The estimate for false positive rate" << endl
        << "   -o, --outDir   STR        Output directory" << endl
        << endl
        << "Optional arguments:" << endl
        << "   -c, --const    FLT        The constant for weights of soft clauses[1000000]" << endl
        // << "   -m, --maxMut   INT        Max number mutations to be eliminated [0]" << endl
        // << "   -b, --bulk     INT        Bulk sequencing file [\"\"]" << endl
        // << "   -e, --delta    FLT        Delta in VAF [0.01]" << endl
        // << "   -v, --truevaf             Use true VAFs instead of noisy VAFs [false]" << endl
        << "   -s, --solver   STR        Name of MaxSAT solver. Available options are:" << endl
        << "                             qmaxsat/maxino/openwbo/aspino/mscg/maxhs [\"openwbo\"]" << endl
        << "                             Note: for another solver, pass the path to binary" << endl
        << "   -i, --integer             Round weights to their nearest integers [false]" << endl
        << "   -z, --coefficient         The coefficient for rounding weights [1000]" << endl
        // << "   -t, --threads  INT        Number of threads [1]" << endl
        << endl
        << "Other arguments:" << endl
        << "   -h, --help                Show this help message and exit" << endl;
}

bool command_line_parser(int argc, char *argv[])
{
    int index;
    char c;
    static struct option longOptions[] = 
    {
    //     {"progress",                no_argument,        &progressRep,       1},
        {"file",                   required_argument,  0,                  'f'},
        {"fnRate",                 required_argument,  0,                  'n'},
        {"fpRate",                 required_argument,  0,                  'p'},
        {"outDir",                 required_argument,  0,                  'o'},
        {"const",                  required_argument,  0,                  'c'},
        // {"maxMut",                 required_argument,  0,                  'm'},
        // {"bulk",                   required_argument,  0,                  'b'},
        // {"delta",                  required_argument,  0,                  'e'},
        // {"vafTrue",                no_argument,        0,                  'v'},
        {"solver",                 required_argument,  0,                  's'},
        {"integer",                required_argument,  0,                  'i'},
        {"coefficient",            required_argument,  0,                  'z'},
        // {"threads",                required_argument,  0,                  't'},
        {"help",                   no_argument,        0,                  'h'},
        {0,0,0,0}
    };

    // while ( (c = getopt_long ( argc, argv, "f:n:p:o:m:b:e:s:t:vh", longOptions, &index))!= -1 )
    while ( (c = getopt_long ( argc, argv, "f:n:p:o:c:s:z:ih", longOptions, &index))!= -1 )
    {
        switch (c)
        {
            case 'f':
                par_inputFile = optarg;
                break;
            case 'n':
                par_fnRate = str2double(optarg);
                if(par_fnRate <= 0 || par_fnRate >= 1)
                {
                    cerr<< "[ERROR] Estimate for false negative rate should be an double in range (0,1)" << endl;
                    return false;
                }
                break;
            case 'p':
                par_fpRate = str2double(optarg);
                if(par_fpRate <= 0 || par_fpRate >= 1)
                {
                    cerr<< "[ERROR] Estimate for false positive rate should be an double in range (0,1)" << endl;
                    return false;
                }
                break;
            case 'o':
                par_outDir = optarg;
                break;
            case 'c':
                par_const = str2int(optarg);
                if(par_const <= 1)
                {
                    cerr<< "[ERROR] Constant for weights should be an integer > 1" << endl;
                    return false;
                }
                break;
            // case 'm':
            //     par_maxColRemove = str2int(optarg);
            //     if(par_maxColRemove < 0)
            //     {
            //         cerr<< "[ERROR] Maximum number of mutation removal should be an integer >= 0" << endl;
            //         return false;
            //     }
            //     break;
            // case 'b':
            //     par_bulkFile = optarg;
            //     break;
            // case 'e':
            //     par_delta = str2double(optarg);
            //     if(par_delta <= 0)
            //     {
            //         cerr<< "[ERROR] Delta should be a floating point number > 0" << endl;
            //         return false;
            //     }
            //     break;
            // case 'v':
            //     par_isTrueVAF = true;
            //     break;
            case 's':
                par_maxSolver = optarg;
                break;
            case 'i':
                INT_WEIGHTS = true;
                break;
            case 'z':
                par_precisionFactor = str2int(optarg);
                if(par_precisionFactor < 1)
                {
                    cerr<< "[ERROR] Rounding coefficient should be an integer >= 1" << endl;
                    return false;
                }
                break;
            // case 't':
            //     par_threads = str2int(optarg);
            //     if(par_threads != 1)
            //     {
            //         cerr<< "[ERROR] Only single thread is supported at the moment!" << endl;
            //         return false;
            //     }
            //     break;
            case 'h':
                print_usage();
                print_help();
                exit(EXIT_SUCCESS);
        }
    }

    if(par_inputFile == "")
    {
        cerr<< "[ERROR] option -f/--file is required" << endl;
        print_usage();
        return false;
    }

    if(par_outDir == "")
    {
        cerr<< "[ERROR] option -o/--outDir is required" << endl;
        print_usage();
        return false;
    }

    if(par_fnRate < 0)
    {
        cerr<< "[ERROR] option -n/--fnRate is required" << endl;
        print_usage();
        return false;
    }

    if(par_fpRate < 0)
    {
        cerr<< "[ERROR] option -p/--fpRate is required" << endl;
        print_usage();
        return false;
    }

    if(INT_WEIGHTS == false)
    {
        par_precisionFactor = 1;
    }

    string exeDir = get_dir_path(get_exe_path());
    if(par_maxSolver == "qmaxsat")
        MAXSAT_EXE = exeDir + "/solvers/qmaxsat/qmaxsat14.04auto-glucose3_static";
    else if(par_maxSolver == "openwbo")
        MAXSAT_EXE = exeDir + "/solvers/open-wbo/open-wbo_glucose4.1_static";
    else if(par_maxSolver == "maxino")
        MAXSAT_EXE = exeDir + "/solvers/maxino/maxino-2015-k16-static";
    else if(par_maxSolver == "aspino")
        MAXSAT_EXE = exeDir + "/solvers/aspino/aspino-static -mode=maxsat";
    else if(par_maxSolver == "mscg")
        MAXSAT_EXE = exeDir + "/solvers/mscg/mscg15b-linux-x86-64";
    else if(par_maxSolver == "maxhs")
        MAXSAT_EXE = exeDir + "/solvers/maxhs/maxhs";
    // else // wrong solver name, use openwbo
    else // expect the full path to the binary of the solver
    {
        MAXSAT_EXE = par_maxSolver;
        // cerr<< "[WARNING] Wrong solver name! Using default solver (openwbo)..." << endl;
        // MAXSAT_EXE = exeDir + "/solvers/open-wbo/open-wbo_glucose4.1_static";
    }

    return true;
}

void get_input_data(string path)
{
    int i, j;
    string tmpStr;
    string line;
    ifstream fin(path.c_str());
    if(fin.is_open() == false)
    {
        cerr<< "Could not open file: " << path << endl;
        exit(EXIT_FAILURE);
    }
    // process the header
    getline(fin, line);
    istringstream sin1(line);
    while(sin1 >> tmpStr)
    {
        mutId.push_back(tmpStr);
    }
    numMut = mutId.size() - 1;
    //
    i = 0;
    while(getline(fin, line))
    {
        istringstream sin(line);
        sin >> tmpStr; // cell name
        cellId.push_back(tmpStr);
        for(int j = 0; j < numMut; j++)
        {
            sin >> mat[i][j];
        }
        i++;
    }
    numCell = i;
    fin.close();
    // // artificial cell and mutation
    // mutId.push_back("mutX");
    // cellId.push_back("cellX");
    // for(j = 0; j < numMut; j++)
    // {
    //     mat[numCell][j] = 0;
    // }
    // for(i = 0; i <= numCell; i++)
    // {
    //     mat[i][numMut] = 1;
    // }
}

void set_y_variables()
{
    int i, j;
    numVarY = 0;

    // for(i = 0; i <= numCell; i++)
    for(i = 0; i < numCell; i++)
    {
        // for(j = 0; j <= numMut; j++)
        for(j = 0; j < numMut; j++)
        {
            numVarY++;
            var_y[i][j] = numVarY;
            map_y2ij[numVarY] = make_pair<int, int>(i, j);
        }
    }
}

void set_x_variables()
{
    int i, j;
    numVarX = 0;

    // for(i = 0; i <= numCell; i++)
    for(i = 0; i < numCell; i++)
    {
        // for(j = 0; j <= numMut; j++)
        for(j = 0; j < numMut; j++)
        {
            numVarX++;
            var_x[i][j] = startVarX + numVarX;
        }
    }
}

void set_b_variables()
{
    int i, j, p, q;
    numVarB = 0;

    // for(p = 0; p <= numMut; p++)
    for(p = 0; p < numMut; p++)
    {
        // for(q = 0; q <= numMut; q++)
        for(q = 0; q < numMut; q++)
        {
            for(i = 0; i < 2; i++)
            {
                for(j = 0; j < 2; j++)
                {
                    numVarB++;
                    var_b[p][q][i][j] = startVarB + numVarB;
                }
            }
        }
    }
}

void set_k_variables()
{
    int p;
    numVarK = 0;

    for(p = 0; p <= numMut; p++)
    {
        numVarK++;
        var_k[p] = startVarK + numVarK;
    }
}

void set_a_variables()
{
    int p, q;
    numVarA = 0;

    for(p = 0; p <= numMut; p++)
    {
        for(q = 0; q <= numMut; q++)
        {
            numVarA++;
            var_a[p][q] = startVarA + numVarA;
            // map_a2pq[startVarA + numVarA] = make_pair<int, int>(p, q);
        }
    }
}

void add_variable_clauses()
{
    int i, j;
    string str_fnWeight;
    string str_fnWeight_neg;
    string str_fpWeight;
    string str_fpWeight_neg;

    numZero = 0;
    numOne = 0;
    numTwo = 0;

    if(INT_WEIGHTS)
    {
        str_fnWeight = int2str(round(par_fnWeight));
        str_fnWeight_neg = int2str(round(par_fnWeight_neg));
        str_fpWeight = int2str(round(par_fpWeight));
        str_fpWeight_neg = int2str(round(par_fpWeight_neg));
    }
    else
    {
        str_fnWeight = double2str(par_fnWeight);
        str_fnWeight_neg = double2str(par_fnWeight_neg);
        str_fpWeight = double2str(par_fpWeight);
        str_fpWeight_neg = double2str(par_fpWeight_neg);
    }

    for(i = 0; i < numCell; i++)
    {
        for(j = 0; j < numMut; j++)
        {
            // fout<< weight_x[map_x2ij[i].first][map_x2ij[i].second] << " " << -1*i << " 0\n";
            if(mat[i][j] == 0)
            {
                numZero++;
                clauseSoft.push_back(str_fnWeight + " " + int2str(var_x[i][j]));
                clauseSoft.push_back(str_fnWeight_neg + " " + int2str(-1*var_x[i][j]));
                clauseHard.push_back(int2str(-1*var_x[i][j]) + " " + int2str(var_y[i][j]));
                clauseHard.push_back(int2str(var_x[i][j]) + " " + int2str(-1*var_y[i][j]));
            }
            else if (mat[i][j] == 1)
            {
                numOne++;
                clauseSoft.push_back(str_fpWeight + " " + int2str(var_x[i][j]));
                clauseSoft.push_back(str_fpWeight_neg + " " + int2str(-1*var_x[i][j]));
                clauseHard.push_back(int2str(var_x[i][j]) + " " + int2str(var_y[i][j]));
                clauseHard.push_back(int2str(-1*var_x[i][j]) + " " + int2str(-1*var_y[i][j]));
            }
            else // mat[i][j] == 2 (not available)
            {
                numTwo++;
                clauseHard.push_back(int2str(-1*var_x[i][j]) + " " + int2str(var_y[i][j]));
                clauseHard.push_back(int2str(var_x[i][j]) + " " + int2str(-1*var_y[i][j]));
            }
        }
    }
}

void add_conflict_clauses()
{
    int i;
    int p, q;
    for(i = 0; i < numCell; i++)
    {
        for(p = 0; p < numMut; p++)
        {
            for(q = p; q < numMut; q++)
            {
                // ~Yip v ~Yiq v Bpq11
                clauseHard.push_back(int2str(-1*var_y[i][p]) + " " + int2str(-1*var_y[i][q]) + " " + int2str(var_b[p][q][1][1]));
                // Yip v ~Yiq v Bpq01
                clauseHard.push_back(int2str(var_y[i][p]) + " " + int2str(-1*var_y[i][q]) + " " + int2str(var_b[p][q][0][1]));
                // ~Yip v Yiq v Bpq10
                clauseHard.push_back(int2str(-1*var_y[i][p]) + " " + int2str(var_y[i][q]) + " " + int2str(var_b[p][q][1][0]));
                // if(par_maxColRemove > 0) // column elimination enabled
                // {
                //     // Kp v Kq v ~Bpq01 v ~Bpq10 v ~Bpq11
                //     clauseHard.push_back(int2str(var_k[p]) + " " + int2str(var_k[q]) + " " + int2str(-1*var_b[p][q][0][1]) + " " + int2str(-1*var_b[p][q][1][0]) + " " + int2str(-1*var_b[p][q][1][1]));
                // }
                // else // column elimination disabled
                // {
                    // ~Bpq01 v ~Bpq10 v ~Bpq11
                    clauseHard.push_back(int2str(-1*var_b[p][q][0][1]) + " " + int2str(-1*var_b[p][q][1][0]) + " " + int2str(-1*var_b[p][q][1][1]));
                // }
            }
        }
    }
}

int next_comb(int comb[], int k, int n)
{
    int i = k - 1;
    ++comb[i];
    while ((i >= 0) && (comb[i] >= n - k + 1 + i))
    {
        --i;
        ++comb[i];
    }

    if (comb[0] > n - k) /* Combination (n-k, n-k+1, ..., n) reached */
        return 0; /* No more combinations can be generated */

    /* comb now looks like (..., x, n, n, n, ..., n).
    Turn it into (..., x, x + 1, x + 2, ...) */
    for (i = i + 1; i < k; ++i)
        comb[i] = comb[i - 1] + 1;

    return 1;
}

void add_column_clauses()
{
    int i;
    // code for C(n, k) 
    // n choose k
    int n = numMut;
    int k = par_maxColRemove + 1;
    int comb[numMut + 10]; // comb[i] is the index of the i-th element in the combination
    for (i = 0; i < k; i++)
        comb[i] = i;

    do
    {
        string tmpClause = "";
        for(i = 0; i < k; i++)
            tmpClause += int2str(-1*var_k[comb[i]]) + " ";
        clauseHard.push_back(tmpClause);
    }while(next_comb(comb, k, n));
}

void add_column_clauses_weight()
{
    int i;
    // int colWeight = numCell / 2;
    // int colWeight = 20;
    string str_colWeight = int2str(par_colWeight);
    for(i = 0; i < numMut; i++)
    {
        clauseSoft.push_back(str_colWeight + " " + int2str(-1*var_k[i]));
    }
}

void add_vaf_clauses()
{
    int t, r;
    int p, q;

    // 1.(a)
    // ~K(numMut)
    if(par_maxColRemove > 0)
    {
        clauseHard.push_back(int2str(-1*var_k[numMut]));   
    }

    // 1.(b)
    // for all rows t, Y(t, numMut) = 1
    for(t = 0; t <= numCell; t++)
    {
        clauseHard.push_back(int2str(var_y[t][numMut]));
    }

    // 1.(c)
    // for all columns p != numMut, Y(numCell, p) = 0
    for(p = 0; p < numMut; p++)
    {
        clauseHard.push_back(int2str(-1*var_y[numCell][p]));
    }

    // // 2.(old): ~a(p,q) v ~a(q,p)
    // for(p = 0; p < numMut; p++)
    // {
    //     for(q = 0; q < numMut; q++)
    //     {
    //         // for all pairs of mutations p and q (including p=q)
    //         clauseHard.push_back(int2str(-1*var_a[p][q]) + " " + int2str(-1*var_a[q][p]));
    //     }
    // }

    // 2.(a)
    // (a(p,q) v a(q,p)) => (~K(p) ^ ~K(q))
    // (~K(p) v ~a(p,q)) ^ (~K(p) v ~a(q,p)) ^ (~a(p,q) v ~K(q)) ^ (~a(q,p) v ~K(q))
    if(par_maxColRemove > 0) // FIXME: should I have this condition or not?
    {
        for(p = 0; p < numMut; p++)
        {
            for(q = 0; q < numMut; q++)
            {
                clauseHard.push_back(int2str(-1*var_k[p]) + " " + int2str(-1*var_a[q][p]));
                clauseHard.push_back(int2str(-1*var_k[p]) + " " + int2str(-1*var_a[q][p]));
                clauseHard.push_back(int2str(-1*var_a[p][q]) + " " + int2str(-1*var_k[q]));
                clauseHard.push_back(int2str(-1*var_a[q][p]) + " " + int2str(-1*var_k[q]));
            }
        }
    }

    // 2.(b)
    // for a given mutation q != x, V_{for all p != q} a(p,q)
    for(q = 0; q < numMut; q++)
    {
        string tmpClause = "";
        for(p = 0; p <= numMut; p++)
        {
            if(p != q)
            {
                tmpClause += int2str(var_a[p][q]) + " ";
            }
        }
        clauseHard.push_back(tmpClause);
    }

    // 2.(c)
    // (a(p,q) ^ Y(t, q)) => (a(p,q) ^ Y(t,p))
    // ~a(p,q) v ~Y(t, q) v Y(t,p)
    for(t = 0; t <= numCell; t++)
    {
        for(p = 0; p <= numMut; p++)
        {
            for(q = 0; q <= numMut; q++)
            {
                clauseHard.push_back(int2str(-1*var_a[p][q]) + " " + int2str(-1*var_y[t][q]) + " " + int2str(var_y[t][p]));
            }
        }
    }

    // 2.(d)
    // a(p,q) => vafP(p,q)
    // ~a(p,q) v vafP(p,q)
    for(p = 0; p <= numMut; p++)
    {
        for(q = 0; q <= numMut; q++)
        {
            if(vafP[p][q] == 0)
            {
                clauseHard.push_back(int2str(-1*var_a[p][q]));
            }
        }
    }

    // 2.(e)
    // cout<< "from " << startVarW << endl;
    numVarW = 0;
    for(p = 0; p <= numMut; p++)
    {
        for(q = 0; q <= numMut; q++)
        {
            // if(p != q) // FIXME: double check
            {
                string tmpClause = "";
                for(t = 0; t <= numCell; t++)
                {
                    numVarW++;
                    clauseHard.push_back(int2str(startVarW+numVarW) + " " + int2str(var_y[t][q]));
                    clauseHard.push_back(int2str(startVarW+numVarW) + " " + int2str(-1*var_y[t][p]));
                    tmpClause += int2str(-1*(startVarW+numVarW)) + " ";
                }
                tmpClause += int2str(var_a[p][q]);
                if(par_maxColRemove > 0)
                {
                    tmpClause += " " + int2str(var_k[p]) + " " + int2str(var_k[q]);   
                }
                clauseHard.push_back(tmpClause);
                // cout<< "new " << tmpClause << endl;
            }
        }
    }

    // 3.
    // (a(p,q) ^ a(p,r) ^ ~a(q,r) ^ ~a(r,q)) => vafT(p,q,r)
    // ~a(p,q) v ~a(p,r) v a(q,r) v a(r,q) v vafT(p,q,r)
    for(p = 0; p <= numMut; p++)
    {
        for(q = 0; q <= numMut; q++)
        {
            for(r = 0; r <= numMut; r++)
            {
                if(q < r) // FIXME: double check
                {
                    if(vafT[p][q][r] == 0)
                    {
                        clauseHard.push_back(int2str(-1*var_a[p][q]) + " " + int2str(-1*var_a[p][r]) + " " + int2str(var_a[q][r]) + " " + int2str(var_a[r][q]));
                    }
                }
            }
        }
    }
}

void write_maxsat_input(string path)
{
    int i, j;

    int64_t hardWeight = numZero * (int64_t)ceil(par_fnWeight) 
                       + numZero * (int64_t)ceil(par_fnWeight_neg) 
                       + numOne  * (int64_t)ceil(par_fpWeight) 
                       + numOne  * (int64_t)ceil(par_fpWeight_neg);
                    //    + numMut * par_colWeight + 1;

    ofstream fout(path.c_str());
    if(fout.is_open() == false)
    {
        cerr<< "Could not open file: " << path << endl;
        exit(EXIT_FAILURE);
    }
    //
    if(IS_PWCNF)
    {
        // fout<< "p wcnf " << numVarY + numVarX + numVarB + numVarK + numVarA + numVarW << " " << clauseSoft.size() + clauseHard.size() << " " << hardWeight << "\n";
        fout<< "p wcnf " << numVarY + numVarX + numVarB << " " << clauseSoft.size() + clauseHard.size() << " " << hardWeight << "\n";
    }
    else
    {
        // fout<< "p wcnf " << numVarY + numVarX + numVarB + numVarK + numVarA + numVarW << " " << clauseSoft.size() + clauseHard.size() << "\n";
        fout<< "p wcnf " << numVarY + numVarX + numVarB << " " << clauseSoft.size() + clauseHard.size() << "\n";
    }
    // soft clauses
    for(i = 0; i < clauseSoft.size(); i++)
    {
        fout<< clauseSoft[i] << " 0\n";
    }
    // hard clauses
    for(i = 0; i < clauseHard.size(); i++)
    {
        fout<< hardWeight << " " << clauseHard[i] << " 0\n";
    }

    fout.close();
}

bool read_maxsat_output_columnElim(string path, int &numRemovedCol, set<int> &removedCol)
{
    numRemovedCol = 0;
    string line;
    bool oLine = false, sLine = false, vLine = false;
    ifstream fin(path.c_str());
    if(fin.is_open() == false)
    {
        cerr<< "Could not open file: " << path << endl;
        exit(EXIT_FAILURE);
    }
    // parse
    while(getline(fin, line))
    {
        if(line[0] == 'o')
        {
            oLine = true;
        }
        if(line[0] == 's')
        {
            sLine = true;
        }
        if(line[0] == 'v')
        {
            vLine = true;
            // update the input matrix
            int tmpVar, tmpVarAbs;
            istringstream sin(line.substr(1));
            while(sin >> tmpVar)
            {
                tmpVarAbs = abs(tmpVar);
                if(startVarK < tmpVarAbs && tmpVarAbs <= startVarK + numVarK) // it is a k variable
                {
                    if(tmpVar > 0) // column to be removed
                    {
                        numRemovedCol++;
                        removedCol.insert(tmpVar - (numVarY + numVarX + numVarB) - 1); // 0-based index
                    }
                }
            }
        }
    }
    fin.close();
    return (oLine && sLine && vLine);
}

bool read_maxsat_output_bitFlips(string path, int &flip, int &flip01, int &flip10, int &flip20, int &flip21, set<int> &removedCol)
{
    flip = 0;
    flip01 = 0;
    flip10 = 0;
    flip20 = 0;
    flip21 = 0;
    string line;
    bool oLine = false, sLine = false, vLine = false;
    ifstream fin(path.c_str());
    if(fin.is_open() == false)
    {
        cerr<< "Could not open file: " << path << endl;
        exit(EXIT_FAILURE);
    }
    // parse
    while(getline(fin, line))
    {
        if(line[0] == 'o')
        {
            oLine = true;
        }
        if(line[0] == 's')
        {
            sLine = true;
        }
        if(line[0] == 'v')
        {
            vLine = true;
            // update the input matrix
            int tmpVar, tmpVarAbs, oldVal;
            istringstream sin(line.substr(1));
            while(sin >> tmpVar)
            {
                tmpVarAbs = abs(tmpVar);
                // if(tmpVarAbs <= numVarY && removedCol.find(tmpVarAbs) == removedCol.end())
                if(tmpVarAbs <= numVarY && removedCol.find(map_y2ij[tmpVarAbs].second) == removedCol.end())
                {
                    oldVal = mat[map_y2ij[tmpVarAbs].first][map_y2ij[tmpVarAbs].second];

                    if(oldVal == 0)
                    {
                        if(tmpVar > 0)
                        {
                            mat[map_y2ij[tmpVarAbs].first][map_y2ij[tmpVarAbs].second] = 1;
                            if(map_y2ij[tmpVarAbs].first != numCell && map_y2ij[tmpVarAbs].second != numMut)
                            {
                                flip++;
                                flip01++;
                            }
                        }
                    }
                    else if(oldVal == 1)
                    {
                        if(tmpVar < 0)
                        {
                            mat[map_y2ij[tmpVarAbs].first][map_y2ij[tmpVarAbs].second] = 0;
                            if(map_y2ij[tmpVarAbs].first != numCell && map_y2ij[tmpVarAbs].second != numMut)
                            {
                                flip++;
                                flip10++;
                            }
                        }
                    }
                    else // oldVal == 2
                    {
                        if(tmpVar < 0)
                        {
                            mat[map_y2ij[tmpVarAbs].first][map_y2ij[tmpVarAbs].second] = 0;
                            if(map_y2ij[tmpVarAbs].first != numCell && map_y2ij[tmpVarAbs].second != numMut)
                            {
                                flip++;
                                flip20++;
                            }
                        }
                        else // tmpVar > 0
                        {
                            mat[map_y2ij[tmpVarAbs].first][map_y2ij[tmpVarAbs].second] = 1;
                            if(map_y2ij[tmpVarAbs].first != numCell && map_y2ij[tmpVarAbs].second != numMut)
                            {
                                flip++;
                                flip21++;
                            }
                        }
                    }
                }
                // // output a variables
                // if(startVarA < tmpVarAbs && tmpVarAbs <= startVarA + numVarA) // it is a A variable
                // {
                //     cout<< "a(" << map_a2pq[tmpVarAbs].first << "," << map_a2pq[tmpVarAbs].second << ") = " << (tmpVar > 0 ? 1 : 0) << endl;
                // }
            }
        }
    }
    fin.close();
    return (oLine && sLine && vLine);
}

void write_output_matrix(string path, set<int> &removedCol)
{
    int i, j;
    ofstream fout(path.c_str());
    // header
    for(i = 0; i < mutId.size(); i++)
    // for(i = 0; i < mutId.size() - 1; i++)
    {
        if(removedCol.find(i-1) == removedCol.end()) // column not removed
            fout<< mutId[i] << "\t";
    }
    fout<< "\n";
    //content
    // for(i = 0; i <= numCell; i++)
    for(i = 0; i < numCell; i++)
    {
        fout<< cellId[i] << "\t";
        // for(j = 0; j <= numMut; j++)
        for(j = 0; j < numMut; j++)
        {
            if(removedCol.find(j) == removedCol.end()) // column not removed
                fout<< mat[i][j] << "\t";
        }
        fout<< "\n";
    }

    fout.close();
}

void get_bulk_data(string path)
{
    int p, q, r;
    string tmpStr;
    double tmpFlt;
    int refCount, mutCount;
    string info;
    string line;
    ifstream fin(path.c_str());
    if(fin.is_open() == false)
    {
        cerr<< "Could not open file: " << path << endl;
        exit(EXIT_FAILURE);
    }
    // get header line
    getline(fin, line);
    // get VAF information
    p = 0;
    while(getline(fin, line))
    {
        if(par_isTrueVAF == false)
        {
            istringstream sin(line);
            sin >> tmpStr;
            sin >> tmpStr;
            sin >> tmpStr;
            sin >> mutCount;
            sin >> refCount;
            vaf[p++] = (double)mutCount/(mutCount + refCount);
        }
        else
        {
            int pos = line.find("trueVAF=");
            vaf[p++] = str2double(line.substr(pos+8, line.find(';', pos) - (pos + 8)));
        }
    }
    // artificial mutation; its VAF should be set to 1
    vaf[numMut] = 5;
    // calc vafP
    for(p = 0; p <= numMut; p++)
    {
        for(q = 0; q <= numMut; q++)
        {
            if(vaf[p]*(1+par_delta) >= vaf[q])
            {
                vafP[p][q] = 1;
            }
            else
            {
                vafP[p][q] = 0;
            }
        }
    }
    // calc vafT
    for(p = 0; p <= numMut; p++)
    {
        for(q = 0; q <= numMut; q++)
        {
            for(r = 0; r <= numMut; r++)
            {
                if(vaf[p]*(1+par_delta) >= vaf[q] + vaf[r])
                {
                    vafT[p][q][r] = 1;
                }
                else
                {
                    vafT[p][q][r] = 0;
                }   
            }
        }
    }
}

bool is_conflict_free()
{
    for(int p = 0; p < numMut; p++)
    {
        for(int q = p + 1; q < numMut; q++)
        {
            bool seen11 = false;
            bool seen01 = false;
            bool seen10 = false;
            for(int r = 0; r < numCell; r++)
            {
                if(mat[r][p] == 1 && mat[r][q] == 1) seen11 = true;
                if(mat[r][p] == 0 && mat[r][q] == 1) seen01 = true;
                if(mat[r][p] == 1 && mat[r][q] == 0) seen10 = true;
            }
            if(seen11 && seen01 && seen10) return false;
        }
    }
    return true;
}

int main(int argc, char *argv[])
{
    if(argc <= 1)
    {
        print_usage();
        exit(EXIT_FAILURE);
    }

    if(command_line_parser(argc, argv) == false)
    {
        exit(EXIT_FAILURE);
    }

    // // calculate integer weights; log in base 10
    // par_fnWeight = round(-1 * par_fnRate + log10(par_const));
    // par_fnWeight_neg = round(log10(1 - pow(2, -1 * par_fnRate)) + log10(par_const));
    // par_fpWeight = round(-1 * par_fpRate + log10(par_const));
    // par_fpWeight_neg = round(log10(1 - pow(2, -1 * par_fpRate)) + log10(par_const));

    // calculate integer weights; log in base 10
    par_fnWeight = par_precisionFactor * log(par_const * par_fnRate);
    par_fnWeight_neg = par_precisionFactor * log(par_const * (1 - par_fnRate));
    par_fpWeight = par_precisionFactor * log(par_const * par_fpRate);
    par_fpWeight_neg = par_precisionFactor * log(par_const * (1 - par_fpRate));

    cout<< "par_fnWeight\t" << par_fnWeight << endl;
    cout<< "par_fnWeight_neg\t" << par_fnWeight_neg << endl;
    cout<< "par_fpWeight\t" << par_fpWeight << endl;
    cout<< "par_fpWeight_neg\t" << par_fpWeight_neg << endl;

    // return 0;

    // create working directory if does not exist
    // FIXME: use a more portable mkdir... int mkdir(const char *path, mode_t mode);
    string cmd = "mkdir -p " + par_outDir;
    system(cmd.c_str());
    string fileName = par_outDir + "/" + get_file_name(par_inputFile, true);

    // // set weights according to the new formulation
    // if(par_maxColRemove > 0)
    // {
    //     // par_colWeight = 0; // for old column elimination formulation
    //     par_colWeight = par_maxColRemove;   
    // }
    // else
    // {
    //     par_colWeight = 0;
    // }

    // double cpuTime = getCpuTime();
    double realTime = getRealTime();

    get_input_data(par_inputFile);

    // if(par_bulkFile != "")
    // {
    //     get_bulk_data(par_bulkFile);
    // }
    // set variables
    set_y_variables();
    set_x_variables();
    set_b_variables();
    // if(par_maxColRemove > 0) // column elimination enabled
    // {
    //     set_k_variables();
    // }
    // if(par_bulkFile != "")
    // {
    //     set_a_variables();
    // }
    // add clauses
    add_variable_clauses();
    add_conflict_clauses();
    // if(par_maxColRemove > 0) // column elimination enabled
    // {
    //     // add_column_clauses();
    //     add_column_clauses_weight();
    // }
    // if(par_bulkFile != "")
    // {
    //     add_vaf_clauses();
    // }
    //
    write_maxsat_input(fileName + ".maxSAT.in");
    
    // run Max-SAT solver
    double maxsatTime = getRealTime();
    cmd = MAXSAT_EXE + " " + fileName + ".maxSAT.in" + " > " + fileName + ".maxSAT.out";
    system(cmd.c_str());
    maxsatTime = getRealTime() - maxsatTime;

    int numFlip = 0;
    int numFlip01 = 0;
    int numFlip10 = 0;
    int numFlip20 = 0;
    int numFlip21 = 0;
    int numRemovedCol = 0;
    set<int> removedCol;

    // if(par_maxColRemove > 0)
    // {
    //     if(read_maxsat_output_columnElim(fileName + ".maxSAT.out", numRemovedCol, removedCol) == false)
    //     {
    //         cerr<< "[ERROR] Max-SAT solver faild!"<< endl;
    //         exit(EXIT_FAILURE);
    //     }
    // }
    //
    if(read_maxsat_output_bitFlips(fileName + ".maxSAT.out", numFlip, numFlip01, numFlip10, numFlip20, numFlip21, removedCol) == false)
    {
        cerr<< "[ERROR] Max-SAT solver faild!"<< endl;
        exit(EXIT_FAILURE);
    }

    // solution is found, save it!
    // write_output_matrix(fileName + ".output", removedCol);
    write_output_matrix(fileName + ".CFMatrix", removedCol);
    // report the log file
    ofstream fLog((fileName + ".log").c_str());
    if(fLog.is_open() == false)
    {
        cerr<< "Could not open file: " << fileName + ".log" << endl;
        exit(EXIT_FAILURE);
    }
    // fLog.precision(6);
    // fLog<< fixed;
    fLog<< "FILE_NAME: " << get_file_name(par_inputFile) << "\n";
    fLog<< "NUM_CELLS(ROWS): " << numCell << "\n";
    fLog<< "NUM_MUTATIONS(COLUMNS): " << numMut << "\n";
    fLog<< "FN_RATE: " << par_fnRate << "\n";
    fLog<< "FN_WEIGHT: " << (INT_WEIGHTS ? (int)round(par_fnWeight) : par_fnWeight) << "\n";
    fLog<< "FN_WEIGHT_NEG: " << (INT_WEIGHTS ? (int)round(par_fnWeight_neg) : par_fnWeight_neg) << "\n";
    fLog<< "FP_RATE: " << par_fpRate << "\n";
    fLog<< "FP_WEIGHT: " << (INT_WEIGHTS ? (int)round(par_fpWeight) : par_fpWeight) << "\n";
    fLog<< "FP_WEIGHT_NEG: " << (INT_WEIGHTS ? (int)round(par_fpWeight_neg) : par_fpWeight_neg) << "\n";
    fLog<< "NUM_THREADS: " << par_threads << "\n";
    fLog<< "MODEL_SOLVING_TIME_SECONDS: " << maxsatTime << "\n";
    fLog<< "RUNNING_TIME_SECONDS: " << getRealTime() - realTime << "\n";
    fLog<< "IS_CONFLICT_FREE: " << (is_conflict_free() ? "YES" : "NO") << "\n"; // FIXME: write the function
    fLog<< "TOTAL_FLIPS_REPORTED: " << numFlip01 + numFlip10 << "\n";
    fLog<< "0_1_FLIPS_REPORTED: " << numFlip01 << "\n";
    fLog<< "1_0_FLIPS_REPORTED: " << numFlip10 << "\n";
    fLog<< "2_0_FLIPS_REPORTED: " << numFlip20 << "\n";
    fLog<< "2_1_FLIPS_REPORTED: " << numFlip21 << "\n";
    fLog<< "MUTATIONS_REMOVED_UPPER_BOUND: " << par_maxColRemove << "\n";
    fLog<< "MUTATIONS_REMOVED_NUM: " << numRemovedCol << "\n";
    fLog<< "MUTATIONS_REMOVED_INDEX: ";
    int ii;
    set<int>::iterator it;
    for(ii = 1, it = removedCol.begin(); it != removedCol.end(); it++, ii++)
    {
        fLog<< (*it)+1 << (ii < removedCol.size() ? "," : "");
    }
    fLog << "\n";
    fLog<< "MUTATIONS_REMOVED_NAME: ";
    for(ii = 1, it = removedCol.begin(); it != removedCol.end(); it++, ii++)
    {
        fLog<< mutId[(*it)+1] << (ii < removedCol.size() ? "," : "");
    }
    fLog << "\n";

    fLog.close();

    if(remove((fileName + ".maxSAT.in").c_str()) != 0 )
        cerr<< "Could not remove file:" << fileName + ".maxSAT.in" << endl;
    if(remove((fileName + ".maxSAT.out").c_str()) != 0 )
        cerr<< "Could not remove file:" << fileName + ".maxSAT.out" << endl;

    return EXIT_SUCCESS;
}
