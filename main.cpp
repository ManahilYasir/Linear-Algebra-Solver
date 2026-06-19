```cpp
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <cmath>
#include <chrono>
#include <cctype>

using namespace std;

// ======= CONFIG =======
const int MAX_N = 50;
const double EPS = 1e-12;
const double VERIFY_EPS = 1e-6;

// ======= UTILITIES =======
string trim(const string &s) {
    int a = 0;
    while (a < (int)s.size() && isspace((unsigned char)s[a])) ++a;
    int b = (int)s.size() - 1;
    while (b >= a && isspace((unsigned char)s[b])) --b;
    return (a <= b ? s.substr(a, b - a + 1) : "");
}

bool stringToInt(const string &s, int &out) {
    string x = trim(s);
    if (x.empty()) return false;
    int i = 0, sign = 1;
    if (x[0] == '+' || x[0] == '-') {
        if (x[0] == '-') sign = -1;
        i = 1;
        if (x.size() == 1) return false;
    }
    long long v = 0;
    for (; i < (int)x.size(); ++i) {
        if (!isdigit((unsigned char)x[i])) return false;
        v = v * 10 + (x[i] - '0');
        if (v > 2000000000LL) return false;
    }
    out = int(sign * v);
    return true;
}

bool stringToDoubleSimple(const string &s, double &out) {
    string x = trim(s);
    if (x.empty()) return false;
    int i = 0, sign = 1;
    if (x[0] == '+' || x[0] == '-') {
        if (x[0] == '-') sign = -1;
        i = 1;
        if (x.size() == 1) return false;
    }
    bool seenDigit = false, seenDot = false;
    double ip = 0.0, fp = 0.0, div = 1.0;
    for (; i < (int)x.size(); ++i) {
        char c = x[i];
        if (isdigit((unsigned char)c)) {
            seenDigit = true;
            if (!seenDot) ip = ip * 10.0 + (c - '0');
            else { fp = fp * 10.0 + (c - '0'); div *= 10.0; }
        } else if (c == '.') {
            if (seenDot) return false;
            seenDot = true;
        } else return false;
    }
    if (!seenDigit) return false;
    out = sign * (ip + fp / div);
    return true;
}

bool stringToDoubleFraction(const string &s, double &out) {
    string x = trim(s);
    if (stringToDoubleSimple(x, out)) return true;
    size_t slash = x.find('/');
    if (slash == string::npos) return false;
    int N, D;
    if (!stringToInt(trim(x.substr(0, slash)), N)) return false;
    if (!stringToInt(trim(x.substr(slash + 1)), D)) return false;
    if (D == 0) return false;
    out = double(N) / double(D);
    return true;
}

double readDouble() {
    string line;
    double v;
    while (true) {
        getline(cin, line);
        if (stringToDoubleFraction(line, v)) return v;
        cout << "Invalid number! Allowed: integer, decimal or fraction a/b. Try again: ";
    }
}

int readIntInRange(int L, int R) {
    string line;
    int v;
    while (true) {
        getline(cin, line);
        if (!stringToInt(line, v)) { cout << "Invalid integer. Enter between " << L << " and " << R << ": "; continue; }
        if (v < L || v > R) { cout << "Out of range. Enter between " << L << " and " << R << ": "; continue; }
        return v;
    }
}

void printClean(double x) {
    if (fabs(x - round(x)) < 1e-9) cout << (int)round(x);
    else cout << fixed << setprecision(6) << x;
}

// ======= DATA STRUCTURES =======
struct StepLog {
    vector<string> steps;
    void push(const string &s) { steps.push_back(s); }
    void show() const { for (auto &s : steps) cout << s << "\n"; }
};

struct Analytics {
    long long rowOps = 0;
    long long rowSwaps = 0;
    double execMs = 0.0;
};

struct SystemData {
    string name;
    int m, n;                       // m equations, n variables
    vector<vector<double>> A;       // size m x n
    vector<double> b;               // size m
};

//  MATRIX HELPERS 
double determinantElim(vector<vector<double>> M) {
    int n = M.size();
    if (n == 0) return 0.0;
    double det = 1.0;
    for (int k = 0; k < n; ++k) {
        int p = k;
        for (int i = k + 1; i < n; ++i)
            if (fabs(M[i][k]) > fabs(M[p][k])) p = i;
        if (fabs(M[p][k]) < EPS) return 0.0;
        if (p != k) {
            swap(M[p], M[k]);
            det = -det;
        }
        det *= M[k][k];
        for (int i = k + 1; i < n; ++i) {
            double f = M[i][k] / M[k][k];
            for (int j = k; j < n; ++j)
                M[i][j] -= f * M[k][j];
        }
    }
    return det;
}

int rankOfMatrix(vector<vector<double>> M) {
    int rows = M.size();
    int cols = M.empty() ? 0 : M[0].size();
    int r = 0;
    for (int c = 0; c < cols && r < rows; ++c) {
        int sel = r;
        for (int i = r; i < rows; ++i)
            if (fabs(M[i][c]) > fabs(M[sel][c])) sel = i;
        if (fabs(M[sel][c]) < EPS) continue;
        swap(M[r], M[sel]);
        for (int i = r + 1; i < rows; ++i) {
            double f = M[i][c] / M[r][c];
            for (int j = c; j < cols; ++j) M[i][j] -= f * M[r][j];
        }
        ++r;
    }
    return r;
}

// SOLVERS (Gaussian & Cramer & Gauss-Jordan) 
bool solveGaussian(const SystemData &sys, vector<double> &x, StepLog &steps, Analytics &ana, string &solutionType) {
    int m = sys.m, n = sys.n;
    // Work with augmented m x (n+1) matrix
    vector<vector<double>> M(m, vector<double>(n + 1));
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) M[i][j] = sys.A[i][j];
        M[i][n] = sys.b[i];
    }

    int r = 0;
    vector<int> pivotCol(m, -1);
    for (int c = 0; c < n && r < m; ++c) {
        int p = r;
        for (int i = r; i < m; ++i)
            if (fabs(M[i][c]) > fabs(M[p][c])) p = i;
        if (fabs(M[p][c]) < EPS) continue;
        if (p != r) {
            swap(M[p], M[r]);
            ana.rowSwaps++;
            steps.push("Swapped row " + to_string(r+1) + " with row " + to_string(p+1));
        }
        pivotCol[r] = c;
        for (int i = r + 1; i < m; ++i) {
            if (fabs(M[i][c]) < EPS) continue;
            double f = M[i][c] / M[r][c];
            for (int j = c; j <= n; ++j) {
                M[i][j] -= f * M[r][j];
                ana.rowOps++;
            }
            steps.push("Eliminated row " + to_string(i+1));
        }
        ++r;
    }

    // Check inconsistency
    for (int i = r; i < m; ++i) {
        double maxA = 0.0;
        for (int j = 0; j < n; ++j) maxA = max(maxA, fabs(M[i][j]));
        if (maxA < EPS && fabs(M[i][n]) > EPS) {
            solutionType = "No Solution";
            return false;
        }
    }

    if (r < n) solutionType = "Infinite Solutions";
    else solutionType = "Unique Solution";

    // Back substitution for one particular solution (free vars=0)
    x.assign(n, 0.0);
    for (int i = r - 1; i >= 0; --i) {
        int pc = pivotCol[i];
        double rhs = M[i][n];
        for (int j = pc + 1; j < n; ++j) rhs -= M[i][j] * x[j];
        if (fabs(M[i][pc]) < EPS) {
            // Shouldn't happen if ranks computed earlier
            x[pc] = 0;
        } else x[pc] = rhs / M[i][pc];
        steps.push("Back-substitution for x[" + to_string(pc) + "]");
    }
    return true;
}

bool solveCramer(const SystemData &sys, vector<double> &x, StepLog &steps, Analytics &ana, string &solutionType) {
    if (sys.m != sys.n) { solutionType = "Cramer's requires square matrix"; return false; }
    int n = sys.n;
    double detA = determinantElim(sys.A);
    if (fabs(detA) < EPS) { solutionType = "Singular matrix (Cramer not applicable)"; return false; }
    x.assign(n, 0.0);
    for (int c = 0; c < n; ++c) {
        vector<vector<double>> M = sys.A;
        for (int i = 0; i < n; ++i) M[i][c] = sys.b[i];
        double detC = determinantElim(M);
        x[c] = detC / detA;
        steps.push("Cramer computed x[" + to_string(c) + "]");
    }
    solutionType = "Unique Solution";
    return true;
}

bool solveGaussJordan(const SystemData &sys, vector<double> &x, StepLog &steps, Analytics &ana, string &solutionType) {
    if (sys.m != sys.n) { solutionType = "Gauss-Jordan requires square matrix"; return false; }
    int n = sys.n;
    vector<vector<double>> M(n, vector<double>(2*n));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) M[i][j] = sys.A[i][j];
        for (int j = 0; j < n; ++j) M[i][n+j] = (i == j) ? 1.0 : 0.0;
    }

    for (int i = 0; i < n; ++i) {
        int p = i;
        for (int r = i; r < n; ++r)
            if (fabs(M[r][i]) > fabs(M[p][i])) p = r;
        if (fabs(M[p][i]) < EPS) { solutionType = "Singular matrix (no inverse)"; return false; }
        if (p != i) { swap(M[p], M[i]); ana.rowSwaps++; steps.push("Swapped row " + to_string(i+1)); }
        double div = M[i][i];
        for (int j = 0; j < 2*n; ++j) { M[i][j] /= div; ana.rowOps++; }
        steps.push("Normalized row " + to_string(i+1));
        for (int r = 0; r < n; ++r) {
            if (r == i) continue;
            double f = M[r][i];
            for (int j = 0; j < 2*n; ++j) { M[r][j] -= f * M[i][j]; ana.rowOps++; }
            steps.push("Eliminated row " + to_string(r+1));
        }
    }

    x.assign(n, 0.0);
    for (int i = 0; i < n; ++i) {
        double s = 0.0;
        for (int j = 0; j < n; ++j) s += M[i][n+j] * sys.b[j];
        x[i] = s;
    }
    solutionType = "Unique Solution";
    return true;
}

//DISPLAY (compact augmented block)
void printAugmentedCompact(const SystemData &sys) {
    cout << "Augmented Matrix [A | b]:\n";
    for (int i = 0; i < sys.m; ++i) {
        for (int j = 0; j < sys.n; ++j) {
            if (fabs(sys.A[i][j] - round(sys.A[i][j])) < 1e-9) 
                cout << setw(8) << (int)round(sys.A[i][j]);
            else 
                cout << setw(8) << fixed << setprecision(6) << sys.A[i][j];
        }
        cout << " | ";
        if (fabs(sys.b[i] - round(sys.b[i])) < 1e-9) 
            cout << setw(8) << (int)round(sys.b[i]);
        else 
            cout << setw(8) << fixed << setprecision(6) << sys.b[i];
        cout << "\n";
    }
}


void printMenuHeader() {
    cout << "\n-----------------------------------\n";
    cout << "-----    LINEAR ALGEBRA SOLVER    ----\n";
    cout << "-----------------------------------\n";
    cout << "1. Add system\n";
    cout << "2. Display systems\n";
    cout << "3. Solve a system\n";
    cout << "4. Update equation\n";
    cout << "5. Delete equation\n";
    cout << "6. Exit\n";
    cout << "Choice: ";
}

void suggestMethods(const SystemData &sys) {
    cout << "\n--- Method suggestions based on system ---\n";
    cout << "- Gaussian elimination: works for any (rectangular or square)\n";
    if (sys.m == sys.n) {
        double detA = determinantElim(sys.A);
        if (fabs(detA) > EPS) {
            cout << "- Cramer's rule: applicable (square & non-singular)\n";
            cout << "- Gauss-Jordan (matrix inversion): applicable\n";
            cout << "Best recommended: Gaussian elimination (numerically stable & general)\n";
        } else {
            cout << "- Matrix appears singular (det ≈ 0). Use Gaussian elimination to detect infinite/no solution\n";
            cout << "Best recommended: Gaussian elimination with pivoting\n";
        }
    } else {
        cout << "- System is rectangular: Gaussian elimination is the appropriate choice\n";
    }
}

//  MAIN 
int main() {
    vector<SystemData> systems;

    while (true) {
        printMenuHeader();
        int choice = readIntInRange(1, 6);
        if (choice == 6) { cout << "Goodbye.\n"; break; }

        if (choice == 1) {                 // Add system
            SystemData sys;
            cout << "Enter system name: ";
            getline(cin, sys.name);

            cout << "1) Square system (n x n)\n2) Rectangular system (m x n)\nChoice: ";
            int t = readIntInRange(1, 2);
            if (t == 1) {
                cout << "Enter n (1-" << MAX_N << "): ";
                sys.n = sys.m = readIntInRange(1, MAX_N);
            } else {
                cout << "Enter number of equations m (1-" << MAX_N << "): ";
                sys.m = readIntInRange(1, MAX_N);
                cout << "Enter number of variables n (1-" << MAX_N << "): ";
                sys.n = readIntInRange(1, MAX_N);
            }

            sys.A.assign(sys.m, vector<double>(sys.n, 0.0));
            sys.b.assign(sys.m, 0.0);

            cout << "Enter matrix A coefficients (row-wise). You may enter decimals or fractions (a/b):\n";
            for (int i = 0; i < sys.m; ++i) {
                for (int j = 0; j < sys.n; ++j) {
                    cout << "A[" << i << "][" << j << "]: ";
                    sys.A[i][j] = readDouble();
                }
            }

            cout << "Enter RHS vector b (one value per equation):\n";
            for (int i = 0; i < sys.m; ++i) {
                cout << "b[" << i << "]: ";
                sys.b[i] = readDouble();
            }

            systems.push_back(sys);
            cout << "System added.\n";
        }

        else if (choice == 2) {           // Display systems
            if (systems.empty()) {
                cout << "No systems stored.\n";
                continue;
            }
            cout << "\nStored systems:\n";
            for (int i = 0; i < (int)systems.size(); ++i) {
                cout << i+1 << ". " << systems[i].name << " (m=" << systems[i].m << ", n=" << systems[i].n << ")\n";
            }
            cout << "Enter index to view (0=back): ";
            int idx = readIntInRange(0, (int)systems.size());
            if (idx == 0) continue;
            printAugmentedCompact(systems[idx-1]);
        }

        else if (choice == 3) {           // Solve a system
            if (systems.empty()) {
                cout << "No systems available.\n";
                continue;
            }
            cout << "Select system index to solve: ";
            int idx = readIntInRange(1, (int)systems.size()) - 1;
            SystemData &sys = systems[idx];

            // show augmented compact
            cout << "System selected: " << sys.name << "\n";
            printAugmentedCompact(sys);

            // compute ranks
            vector<vector<double>> Acopy = sys.A;
            int rA = rankOfMatrix(Acopy);
            vector<vector<double>> aug = sys.A;
            for (int i = 0; i < sys.m; ++i) aug[i].push_back(sys.b[i]);
            int rAug = rankOfMatrix(aug);

            cout << "Rank(A) = " << rA << ", Rank(A|b) = " << rAug << "\n";

            suggestMethods(sys);

            // New prompt for recommended vs all methods
            cout << "\nHow do you want to solve this system?\n";
            cout << "1) Solve using the recommended method\n";
            cout << "2) Solve using all applicable methods\nChoice: ";
            int mode = readIntInRange(1, 2);

            cout << "Show steps?\n1) Answer only\n2) Step-by-step\nChoice: ";
            int showSteps = readIntInRange(1, 2);

            // Build applicable methods list
            vector<int> applicable;
            applicable.push_back(1); // Gaussian always possible
            if (sys.m == sys.n) {
                double detA = determinantElim(sys.A);
                if (fabs(detA) > EPS) {
                    applicable.push_back(2); // Cramer
                    applicable.push_back(3); // Gauss-Jordan
                }
            }

            // If user chose recommended method, only keep Gaussian
            if (mode == 1) {
                applicable.clear();
                applicable.push_back(1);
            }

            // Execute selected methods
            for (int mth : applicable) {
                vector<double> x;
                StepLog slog;
                Analytics ana;
                string solType;

                auto t0 = chrono::high_resolution_clock::now();
                bool ok = false;
                if (mth == 1) ok = solveGaussian(sys, x, slog, ana, solType);
                else if (mth == 2) ok = solveCramer(sys, x, slog, ana, solType);
                else if (mth == 3) ok = solveGaussJordan(sys, x, slog, ana, solType);
                auto t1 = chrono::high_resolution_clock::now();
                ana.execMs = chrono::duration<double, milli>(t1 - t0).count();

                cout << "\n----- Method: " << (mth==1?"Gaussian":mth==2?"Cramer's":"Gauss-Jordan") << " -----\n";
                cout << "Solution type: " << solType << "\n";
                if (ok) {
                    cout << "Solution (one particular if infinite):\n";
                    for (int i = 0; i < sys.n; ++i) {
                        cout << "x[" << i << "] = ";
                        printClean(x[i]);
                        cout << "\n";
                    }
                } else {
                    cout << "Method failed or not applicable for this system.\n";
                }

                // analytics
                cout << "--- Analytics ---\n";
                if (sys.m == sys.n) {
                    double det = determinantElim(sys.A);
                    cout << "Determinant: ";
                    if (fabs(det) < EPS) cout << "0 (or near zero)\n";
                    else { printClean(det); cout << "\n"; }
                }
                cout << "Rank(A) = " << rA << "\n";
                cout << "Rank(A|b) = " << rAug << "\n";
                cout << "Row ops (approx): " << ana.rowOps << "\n";
                cout << "Row swaps: " << ana.rowSwaps << "\n";
                cout << "Time (ms): " << ana.execMs << "\n";

                if (showSteps == 2) {
                    cout << "\n--- Steps ---\n";
                    slog.show();
                }

                // verification
                if (ok) {
                    bool v = true;
                    if ((int)x.size() == sys.n) {
                        for (int i = 0; i < sys.m; ++i) {
                            double s = 0.0;
                            for (int j = 0; j < sys.n; ++j) 
                                s += sys.A[i][j] * x[j];
                            if (fabs(s - sys.b[i]) > VERIFY_EPS) { v = false; break; }
                        }
                    } else v = false;

                    if (v) cout << "AX = b  OK\n";
                    else cout << "AX != b  FAILED\n";
                }
            }
        }

        else if (choice == 4) {          // Update equation
            if (systems.empty()) { cout << "No systems.\n"; continue; }
            cout << "Select system index: ";
            int idx = readIntInRange(1, (int)systems.size()) - 1;
            SystemData &sys = systems[idx];
            printAugmentedCompact(sys);

            cout << "Equation index to update (1-" << sys.m << "): ";
            int r = readIntInRange(1, sys.m) - 1;

            cout << "1) Update coefficient in A\n2) Update RHS b\nChoice: ";
            int op = readIntInRange(1, 2);

            if (op == 1) {
                cout << "Column index (1-" << sys.n << "): ";
                int c = readIntInRange(1, sys.n) - 1;
                cout << "New value for A[" << r << "][" << c << "]: ";
                sys.A[r][c] = readDouble();
                cout << "Updated A[" << r << "][" << c << "].\n";
            } else {
                cout << "New value for b[" << r << "]: ";
                sys.b[r] = readDouble();
                cout << "Updated b[" << r << "].\n";
            }
        }

        else if (choice == 5) {          // Delete equation with undo option
            if (systems.empty()) { cout << "No systems.\n"; continue; }
            cout << "Select system index: ";
            int idx = readIntInRange(1, (int)systems.size()) - 1;
            SystemData &sys = systems[idx];
            printAugmentedCompact(sys);

            cout << "Equation index to delete (1-" << sys.m << "): ";
            int e = readIntInRange(1, sys.m) - 1;
            vector<double> backupRow = sys.A[e];
            double backupB = sys.b[e];

            // perform deletion
            sys.A.erase(sys.A.begin() + e);
            sys.b.erase(sys.b.begin() + e);
            sys.m--;

            cout << "Equation removed. Undo delete? 1=Yes  2=No : ";
            int undo = readIntInRange(1, 2);
            if (undo == 1) {
                sys.A.insert(sys.A.begin() + e, backupRow);
                sys.b.insert(sys.b.begin() + e, backupB);
                sys.m++;
                cout << "Deletion undone.\n";
            } else {
                cout << "Equation permanently deleted.\n";
            }
        }
    }

    return 0;
}

