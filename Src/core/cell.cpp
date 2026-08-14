#include "pch.h"
#include "cell.h"
#include "convenience.h"
#include "nos_math.h"
#include <cctype>

namespace {
    /**
     * @brief Evaluates the numerical factor of a single term of a symmetry operation.
     *
     * Understands an empty string (the implicit 1 of "x"), a plain number ("2", "0.5")
     * and a fraction ("1/2").
     *
     * @param number the factor as written in the CIF, without its sign
     * @param value [out] the evaluated factor
     * @return whether the string could be evaluated
     */
    bool eval_symop_factor(const std::string &number, double &value) {
        if (number.empty()) {
            value = 1.0;
            return true;
        }
        auto to_double = [](const std::string &s, double &out) {
            if (s.empty())
                return false;
            try {
                size_t used = 0;
                out = std::stod(s, &used);
                return used == s.length();
            }
            catch (const std::exception &) {
                return false;
            }
        };
        const size_t slash = number.find('/');
        if (slash == std::string::npos)
            return to_double(number, value);
        double numerator = 0.0, denominator = 0.0;
        if (!to_double(number.substr(0, slash), numerator))
            return false;
        if (!to_double(number.substr(slash + 1), denominator))
            return false;
        if (denominator == 0.0)
            return false;
        value = numerator / denominator;
        return true;
    }
}

void cell::parse_symop(const std::string &operation,
                       const std::filesystem::path &filename,
                       int rot[3][3],
                       double translation[3],
                       std::ostream &file) {
    const std::string where = " of symmetry operation \"" + operation + "\" in " + filename.string() + "!";
    // Split into the three comma separated components, dropping any whitespace
    svec components(3);
    int column = 0;
    for (const char c : operation) {
        if (c == ',') {
            column++;
            err_checkf(column < 3, "Found more than 3 comma separated components" + where, file);
        }
        else if (!std::isspace(static_cast<unsigned char>(c)))
            components[column].push_back(c);
    }
    err_checkf(column == 2, "Expected 3 comma separated components" + where, file);

    for (int comp = 0; comp < 3; comp++) {
        const std::string &s = components[comp];
        err_checkf(!s.empty(), "Component " + std::to_string(comp + 1) + " is empty" + where, file);
        rot[comp][0] = rot[comp][1] = rot[comp][2] = 0;
        translation[comp] = 0.0;
        // Walk the signed terms, each of which is either an axis (optionally scaled) or a translation
        size_t pos = 0;
        while (pos < s.length()) {
            double sign = 1.0;
            if (s[pos] == '+')
                pos++;
            else if (s[pos] == '-')
                sign = -1.0, pos++;
            size_t end = pos;
            while (end < s.length() && s[end] != '+' && s[end] != '-')
                end++;
            err_checkf(end > pos, "Found an empty term" + where, file);
            std::string term = s.substr(pos, end - pos);
            pos = end;
            // Pull out the axis name, if this term has one; what remains is its factor
            int axis = -1;
            for (size_t k = 0; k < term.length() && axis == -1; k++) {
                switch (term[k]) {
                case 'x': case 'X': axis = 0; break;
                case 'y': case 'Y': axis = 1; break;
                case 'z': case 'Z': axis = 2; break;
                default: continue;
                }
                term.erase(k, 1);
            }
            // "2*x" carries the same information as "2x"
            term.erase(std::remove(term.begin(), term.end(), '*'), term.end());
            double factor = 0.0;
            err_checkf(eval_symop_factor(term, factor), "Could not interpret the factor \"" + term + "\"" + where, file);
            if (axis == -1) {
                translation[comp] += sign * factor;
                continue;
            }
            const double coefficient = sign * factor;
            const int rounded = static_cast<int>(std::lround(coefficient));
            err_checkf(std::abs(coefficient - rounded) < 1e-6,
                       "The rotation coefficient " + std::to_string(coefficient) + " is not an integer" + where, file);
            rot[comp][axis] += rounded;
        }
    }
}

void cell::eval_symm(std::vector<asym_atom> &asym_atoms) {
    const int ncen = asym_atoms.size();
    vec pos(3);
    vec new_pos(3);
    auto sym_ = sym;
    auto trans_ = trans;
    int num_sym = trans[0].size();
    int idx = 0;
    for (asym_atom a : asym_atoms) {
        int count = 0;
        pos[0] = a.frac_pos[0];
        pos[1] = a.frac_pos[1];
        pos[2] = a.frac_pos[2];
        for (int t = 0; t < num_sym; t++) {
            const vec trans_temp = { trans[0][t], trans[1][t], trans[2][t] };
            const vec2 rot_temp = { { (double)sym[0][0][t], (double)sym[0][1][t], (double)sym[0][2][t] },
                                     { (double)sym[1][0][t], (double)sym[1][1][t], (double)sym[1][2][t] },
                                     { (double)sym[2][0][t], (double)sym[2][1][t], (double)sym[2][2][t] } };
            vec temp_pos = self_dot(rot_temp, pos, false);
            new_pos[0] = temp_pos[0] + trans_temp[0];
            new_pos[1] = temp_pos[1] + trans_temp[1];
            new_pos[2] = temp_pos[2] + trans_temp[2];
            if (check_special(pos, new_pos)) {
                count++;
            }
        }
        asym_atoms[idx].asym_fact = 1.0 / count;
        idx++;
    }
    //closing function
}

bool cell::check_special(const vec &pos1, const vec &pos2) {
    bool special = true;
    double dist1 = (pos2[0] - pos1[0]);
    double dist2 = (pos2[1] - pos1[1]);
    double dist3 = (pos2[2] - pos1[2]);
    double mod1 = std::fmod(dist1, 1);
    double mod2 = std::fmod(dist2, 1);
    double mod3 = std::fmod(dist3, 1);
    if (std::abs(mod1) > 1e-10 || std::abs(mod2) > 1e-10 || std::abs(mod3) > 1e-10) {
        special = false;
    }
    return special;
    //closing function
}