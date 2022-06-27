#define _CRT_NONSTDC_NO_WARNINGS
#include <bits/stdc++.h>
#include <random>
#include <unordered_set>
#ifdef _MSC_VER
#define ENABLE_VIS
#define ENABLE_DUMP
#include <conio.h>
#include <ppl.h>
#include <filesystem>
#ifdef ENABLE_VIS
#include <opencv2/core.hpp>
#include <opencv2/core/utils/logger.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/highgui.hpp>
#endif
#endif

/** compro_io **/

/* tuple */
// out
namespace aux {
    template<typename T, unsigned N, unsigned L>
    struct tp {
        static void output(std::ostream& os, const T& v) {
            os << std::get<N>(v) << ", ";
            tp<T, N + 1, L>::output(os, v);
        }
    };
    template<typename T, unsigned N>
    struct tp<T, N, N> {
        static void output(std::ostream& os, const T& v) { os << std::get<N>(v); }
    };
}
template<typename... Ts>
std::ostream& operator<<(std::ostream& os, const std::tuple<Ts...>& t) {
    os << '[';
    aux::tp<std::tuple<Ts...>, 0, sizeof...(Ts) - 1>::output(os, t);
    return os << ']';
}

template<class Ch, class Tr, class Container>
std::basic_ostream<Ch, Tr>& operator<<(std::basic_ostream<Ch, Tr>& os, const Container& x);

/* pair */
// out
template<class S, class T>
std::ostream& operator<<(std::ostream& os, const std::pair<S, T>& p) {
    return os << "[" << p.first << ", " << p.second << "]";
}
// in
template<class S, class T>
std::istream& operator>>(std::istream& is, std::pair<S, T>& p) {
    return is >> p.first >> p.second;
}

/* container */
// out
template<class Ch, class Tr, class Container>
std::basic_ostream<Ch, Tr>& operator<<(std::basic_ostream<Ch, Tr>& os, const Container& x) {
    bool f = true;
    os << "[";
    for (auto& y : x) {
        os << (f ? "" : ", ") << y;
        f = false;
    }
    return os << "]";
}
// in
template <
    class T,
    class = decltype(std::begin(std::declval<T&>())),
    class = typename std::enable_if<!std::is_same<T, std::string>::value>::type
>
    std::istream& operator>>(std::istream& is, T& a) {
    for (auto& x : a) is >> x;
    return is;
}

/* struct */
template<typename T>
auto operator<<(std::ostream& out, const T& t) -> decltype(out << t.stringify()) {
    out << t.stringify();
    return out;
}

/* setup */
struct IOSetup {
    IOSetup(bool f) {
        if (f) { std::cin.tie(nullptr); std::ios::sync_with_stdio(false); }
        std::cout << std::fixed << std::setprecision(15);
    }
} iosetup(true);

/** string formatter **/
template<typename... Ts>
std::string format(const std::string& f, Ts... t) {
    size_t l = std::snprintf(nullptr, 0, f.c_str(), t...);
    std::vector<char> b(l + 1);
    std::snprintf(&b[0], l + 1, f.c_str(), t...);
    return std::string(&b[0], &b[0] + l);
}

template<typename T>
std::string stringify(const T& x) {
    std::ostringstream oss;
    oss << x;
    return oss.str();
}

/* dump */
#ifdef ENABLE_DUMP
#define DUMPOUT std::cerr
std::ostringstream DUMPBUF;
#define dump(...) do{DUMPBUF<<"  ";DUMPBUF<<#__VA_ARGS__<<" :[DUMP - "<<__LINE__<<":"<<__FUNCTION__<<"]"<<std::endl;DUMPBUF<<"    ";dump_func(__VA_ARGS__);DUMPOUT<<DUMPBUF.str();DUMPBUF.str("");DUMPBUF.clear();}while(0);
void dump_func() { DUMPBUF << std::endl; }
template <class Head, class... Tail> void dump_func(Head&& head, Tail&&... tail) { DUMPBUF << head; if (sizeof...(Tail) == 0) { DUMPBUF << " "; } else { DUMPBUF << ", "; } dump_func(std::move(tail)...); }
#else
#define dump(...) void(0);
#endif

/* timer */
class Timer {
    double t = 0, paused = 0, tmp;
public:
    Timer() { reset(); }
    static double time() {
#ifdef _MSC_VER
        return __rdtsc() / 2.8e9;
#else
        unsigned long long a, d;
        __asm__ volatile("rdtsc"
            : "=a"(a), "=d"(d));
        return (d << 32 | a) / 2.8e9;
#endif
    }
    void reset() { t = time(); }
    void pause() { tmp = time(); }
    void restart() { paused += time() - tmp; }
    double elapsed_ms() const { return (time() - t - paused) * 1000.0; }
} g_timer;

/* rand */
struct Xorshift {
    static constexpr uint64_t M = INT_MAX;
    static constexpr double e = 1.0 / M;
    uint64_t x = 88172645463325252LL;
    Xorshift() {}
    Xorshift(uint64_t seed) { reseed(seed); }
    inline void reseed(uint64_t seed) { x = 0x498b3bc5 ^ seed; for (int i = 0; i < 20; i++) next(); }
    inline uint64_t next() { x = x ^ (x << 7); return x = x ^ (x >> 9); }
    inline int next_int() { return next() & M; }
    inline int next_int(int mod) { return next() % mod; }
    inline int next_int(int l, int r) { return l + next_int(r - l + 1); }
    double next_double() { return next_int() * e; }
} rnd;

/* shuffle */
template<typename T>
void shuffle_vector(std::vector<T>& v, Xorshift& rnd) {
    int n = v.size();
    for (int i = n - 1; i >= 1; i--) {
        int r = rnd.next_int(i);
        std::swap(v[i], v[r]);
    }
}

/* split */
std::vector<std::string> split(std::string str, const std::string& delim) {
    for (char& c : str) if (delim.find(c) != std::string::npos) c = ' ';
    std::istringstream iss(str);
    std::vector<std::string> parsed;
    std::string buf;
    while (iss >> buf) parsed.push_back(buf);
    return parsed;
}

template<typename A, size_t N, typename T> inline void Fill(A(&array)[N], const T& val) {
    std::fill((T*)array, (T*)(array + N), val);
}

template<typename T, typename ...Args> auto make_vector(T x, int arg, Args ...args) { if constexpr (sizeof...(args) == 0)return std::vector<T>(arg, x); else return std::vector(arg, make_vector<T>(x, args...)); }
template<typename T> bool chmax(T& a, const T& b) { if (a < b) { a = b; return true; } return false; }
template<typename T> bool chmin(T& a, const T& b) { if (a > b) { a = b; return true; } return false; }

using std::vector, std::string;
using std::cin, std::cout, std::cerr, std::endl;
using ll = long long;
using pii = std::pair<int, int>;




#define BATCH_TEST
#ifdef BATCH_TEST
#undef ENABLE_DUMP
#endif

constexpr int di8[] = { 0, -1, -1, -1, 0, 1, 1, 1 };
constexpr int dj8[] = { 1, 1, 0, -1, -1, -1, 0, 1 };
using cell_t = short;
constexpr cell_t WALL = -1;
constexpr cell_t NONE = 0;

struct Input;
using InputPtr = std::shared_ptr<Input>;
struct Input {
    int N;
    vector<string> S;
    Input(std::istream& in) {
        in >> N;
        S.resize(N, string(N, ' '));
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                in >> S[i][j];
            }
        }
    }
    string stringify() const {
        string res;
        res += std::to_string(N) + '\n';
        for (const auto& s : S) res += s + '\n';
        return res;
    }
};



uint64_t g_hash_table[42][42][1600];

struct HashSetup {
    HashSetup() {
        std::mt19937_64 engine;
        for (int i = 0; i < 42; i++) {
            for (int j = 0; j < 42; j++) {
                for (int k = 0; k < 1600; k++) {
                    g_hash_table[i][j][k] = engine();
                }
            }
        }
    }
} hash_setup;

struct State;
using StatePtr = std::shared_ptr<State>;
struct State {
    int N;
    int largest;
    int score;
    uint64_t hash;
    std::array<std::array<cell_t, 42>, 42> grid;
    //vector<vector<cell_t>> grid;
    State(InputPtr input) :
        N(input->N), largest(0), score(0), hash(0) {
        for (int i = 0; i < 42; i++) for (int j = 0; j < 42; j++) grid[i][j] = WALL;
        for (int i = 1; i <= N; i++) {
            for (int j = 1; j <= N; j++) {
                grid[i][j] = input->S[i - 1][j - 1] == '.' ? NONE : WALL;
            }
        }
    }
    int calc_space() const {
        int num_space = 0;
        for (int i = 1; i <= N; i++) {
            for (int j = 1; j <= N; j++) {
                num_space += grid[i][j] == NONE;
            }
        }
        return num_space;
    }
    int can_move(int i, int j) const {
        if (grid[i][j]) return 0;
        int sum = 0, ctr = 0;
        for (int d = 0; d < 8; d++) {
            int ni = i + di8[d], nj = j + dj8[d];
            ctr += grid[ni][nj] > 0;
            sum += std::max(cell_t(0), grid[ni][nj]);
        }
        if (!sum && !ctr) return 1;
        if (ctr <= 1) return 0;
        return sum <= largest + 1 ? sum : 0;
    }
    inline uint64_t calc_hash(int i, int j, int num) const {
        return hash ^ g_hash_table[i][j][num];
    }
    bool do_best_move() {
        int best_num = -1, best_i = -1, best_j = -1;
        for (int i = 1; i <= N; i++) {
            for (int j = 1; j <= N; j++) {
                int num = -1;
                if (num = can_move(i, j)) {
                    if (chmax(best_num, num)) {
                        best_i = i;
                        best_j = j;
                    }
                }
            }
        }
        if (best_num == -1) return false;
        move(best_i, best_j, best_num);
        return true;
    }
    void move(int i, int j, int num) {
        chmax(largest, num);
        score += num;
        grid[i][j] = num;
        hash ^= g_hash_table[i][j][num];
    }
    vector<StatePtr> enum_all_next_states() {
        vector<StatePtr> res;
        for (int i = 1; i <= N; i++) {
            for (int j = 1; j <= N; j++) {
                int num = -1;
                if (num = can_move(i, j)) {
                    StatePtr ns = std::make_shared<State>(*this);
                    ns->move(i, j, num);
                    res.push_back(ns);
                }
            }
        }
        return res;
    }
    string stringify() const {
        string res;
        res += format("largest=%d, score=%d\n", largest, score);
        string line = "+";
        for (int i = 0; i < N; i++) line += "---+";
        line += '\n';
        res += line;
        for (int i = 1; i <= N; i++) {
            res += '|';
            for (int j = 1; j <= N; j++) {
                res += grid[i][j] == WALL ? "###|" : (grid[i][j] == NONE ? "   |" : format("%3d|", grid[i][j]));
            }
            res += '\n';
            res += line;
        }
        return res;
    }
    void output(std::ostream& out) const {
        int placed = 0;
        vector<vector<pii>> buckets(largest + 1);
        for (int i = 1; i <= N; i++) {
            for (int j = 1; j <= N; j++) {
                if (grid[i][j] > 0) {
                    buckets[grid[i][j]].emplace_back(i - 1, j - 1);
                    placed++;
                }
            }
        }
        out << placed << '\n';
        for (int n = 1; n <= largest; n++) {
            for (const auto [r, c] : buckets[n]) {
                out << r << ' ' << c << '\n';
            }
        }
    }
};

std::ostream& operator<<(std::ostream& o, const StatePtr& s) {
    o << s->stringify();
    return o;
}

StatePtr greedy(StatePtr state) {
    while (state->do_best_move());
    return state;
}

StatePtr beam_search(StatePtr init_state, int beam_width = 500) {

    struct Cmp {
        bool operator()(const StatePtr& a, const StatePtr& b) const {
            return a->score == b->score ? a->hash < b->hash : a->score > b->score;
        }
    };

    std::set<StatePtr, Cmp> now_states({ init_state });
    StatePtr best_state = init_state;
    const int N = init_state->N;
    while (true) {
        if (now_states.empty()) break;
        std::set<StatePtr, Cmp> next_states;
        for (int i = 0; i < beam_width; i++) {
            if (now_states.empty()) break;
            auto now_state = *now_states.begin(); now_states.erase(now_states.begin());
            for (int r = 1; r <= N; r++) {
                for (int c = 1; c <= N; c++) {
                    int num = -1;
                    if (num = now_state->can_move(r, c)) {
                        int score = now_state->score + num;
                        if (next_states.size() < beam_width) {
                            StatePtr ns = std::make_shared<State>(*now_state);
                            ns->move(r, c, num);
                            next_states.insert(ns);
                        }
                        else if ((*std::prev(next_states.end()))->score < score) {
                            StatePtr ns = std::make_shared<State>(*now_state);
                            ns->move(r, c, num);
                            next_states.insert(ns);
                        }
                    }
                }
            }
        }
        if (next_states.empty()) break;
        now_states.swap(next_states);
        if (best_state->score < (*now_states.begin())->score) {
            best_state = *now_states.begin();
            dump(best_state->score);
        }
    }

    return best_state;
}



int main(int argc, char** argv) {

    Timer timer;

#ifdef HAVE_OPENCV_HIGHGUI
    cv::utils::logging::setLogLevel(cv::utils::logging::LogLevel::LOG_LEVEL_SILENT);
#endif

#if 0
    std::ifstream ifs(R"(C:\Users\komori3\OneDrive\dev\compro\heuristic\tasks\MM137\tester\in\2.in)");
    std::ofstream ofs(R"(C:\Users\komori3\OneDrive\dev\compro\heuristic\tasks\MM137\tester\out\2.out)");
    std::istream& in = ifs;
    std::ostream& out = ofs;
#else
    std::istream& in = cin;
    std::ostream& out = cout;
#endif

    auto input = std::make_shared<Input>(in);

    dump(*input);

    auto state = std::make_shared<State>(input);
    auto best_state = beam_search(state);
    //auto best_state = greedy(state);

    dump(best_state->score);
    best_state->output(out);

    dump(timer.elapsed_ms());

    return 0;
}