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
    double calc_sparseness() const {
        int space = 0;
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                space += S[i][j] == '.';
            }
        }
        return space / double(N * N);
    }
    string stringify() const {
        string res;
        res += std::to_string(N) + '\n';
        for (const auto& s : S) res += s + '\n';
        return res;
    }
};

uint64_t g_hash_table[42*42][40*40];
int g_beamwidth[] = {
    -1,-1,-1,-1,-1,-1,-1,-1,30000,25000,
    20000,15000,12000,11000,10000,9000,8000,7000,6000,5000,
    4500,4000,3500,3000,2700,2500,2200,1900,1700,1500,
    1400,1300,1200,1100,1000,950,900,850,800,750,
    700
};

struct HashSetup {
    HashSetup() {
        std::mt19937_64 engine;
        for (int p = 0; p < 42*42; p++) {
            for (int t = 0; t < 40*40; t++) {
                g_hash_table[p][t] = engine();
            }
        }
    }
};

template<int N>
struct Solver {

    using cell_t = short;
    using ctr_t = uint8_t;
    static constexpr cell_t WALL = -1;
    static constexpr cell_t NONE = 0;

    static constexpr int MAX_T = N * N;
    static constexpr int WP = N + 2;
    static constexpr int SP = WP * (N + 2);
    static constexpr int d8[] = { 1, 1 - WP, -WP, -1 - WP, -1, -1 + WP, WP, 1 + WP };

    static inline int enc(int i, int j) { return i * WP + j; }
    static inline std::pair<int, int> dec(int p) { return { p / WP, p % WP }; }

    struct State;
    using StatePtr = std::shared_ptr<State>;
    struct State {

        int largest;
        int score;
        uint64_t hash;
        std::array<cell_t, SP> grid;
        std::array<ctr_t, SP> ctrs;
        std::array<cell_t, SP> sums;

        State(InputPtr input) :
            largest(0), score(0), hash(0) {
            std::fill(grid.begin(), grid.end(), WALL);
            std::fill(ctrs.begin(), ctrs.end(), 0);
            std::fill(sums.begin(), sums.end(), 0);
            for (int i = 1; i <= N; i++) {
                for (int j = 1; j <= N; j++) {
                    grid[enc(i, j)] = input->S[i - 1][j - 1] == '.' ? NONE : WALL;
                }
            }
        }

        int can_move(int p) const {
            if (grid[p]) return 0;
            if (!(sums[p] | ctrs[p])) return 1;
            if (ctrs[p] <= 1) return 0;
            return sums[p] <= largest + 1 ? sums[p] : 0;
        }

        void move(int p, int num) {
            chmax(largest, num);
            score += num;
            grid[p] = num;
            for (int d = 0; d < 8; d++) {
                int np = p + d8[d];
                ctrs[np]++;
                sums[np] += num;
            }
            hash ^= g_hash_table[p][num];
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
                    int p = enc(i, j);
                    res += grid[p] == WALL ? "###|" : (grid[p] == NONE ? "   |" : format("%3d|", grid[p]));
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
                    int p = enc(i, j);
                    if (grid[p] > 0) {
                        buckets[grid[p]].emplace_back(i - 1, j - 1);
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

        friend std::ostream& operator<<(std::ostream& o, const StatePtr& s) {
            o << s->stringify();
            return o;
        }

    };

    StatePtr beam_search(StatePtr init_state, int beam_width = 1000, double duration = 9500) {

        Timer timer;

        struct Cmp {
            bool operator()(const StatePtr& a, const StatePtr& b) const {
                return a->score == b->score ? a->hash < b->hash : a->score > b->score;
            }
        };

        std::set<StatePtr, Cmp> now_states({ init_state });
        auto best_state = init_state;
        int turn = 0;
        while (true) {
            if (now_states.empty()) break;
            turn++;
            std::set<StatePtr, Cmp> next_states;
            for (int i = 0; i < beam_width; i++) {
                if (now_states.empty()) break;
                auto now_state = *now_states.begin(); now_states.erase(now_states.begin());
                for (int r = 1; r <= N; r++) {
                    for (int c = 1; c <= N; c++) {
                        int p = enc(r, c), num = -1;
                        if (num = now_state->can_move(p)) {
                            int score = now_state->score + num;
                            if (next_states.size() < beam_width) {
                                auto ns = std::make_shared<State>(*now_state);
                                ns->move(p, num);
                                next_states.insert(ns);
                            }
                            else if ((*std::prev(next_states.end()))->score < score) {
                                auto ns = std::make_shared<State>(*now_state);
                                ns->move(p, num);
                                next_states.insert(ns);
                            }
                        }
                    }
                }
            }
            if (next_states.empty() || timer.elapsed_ms() > duration) break;
            now_states.swap(next_states);
            if (best_state->score < (*now_states.begin())->score) {
                best_state = *now_states.begin();
                dump(turn, timer.elapsed_ms(), best_state->score);
            }
        }

        return best_state;
    }

    Solver(InputPtr input, std::ostream& out) {

        Timer timer;

        HashSetup();

        auto state = std::make_shared<State>(input);

        const double duration = 9500 - timer.elapsed_ms();
        int beam_width = (int)round(g_beamwidth[N] / input->calc_sparseness());
        auto best_state = beam_search(state, beam_width, duration);

        dump(best_state->score);
        best_state->output(out);

        dump(timer.elapsed_ms());

    }

};

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

#define X(M) { if (input->N == M) Solver<M>(input, out); }
    X(8); X(9); X(10);
    X(11); X(12); X(13); X(14); X(15); X(16); X(17); X(18); X(19); X(20);
    X(21); X(22); X(23); X(24); X(25); X(26); X(27); X(28); X(29); X(30);
    X(31); X(32); X(33); X(34); X(35); X(36); X(37); X(38); X(39); X(40);
#undef X

    return 0;
}