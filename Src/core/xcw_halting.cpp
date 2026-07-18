#include "pch.h"
#include "xcw_halting.h"

double std_normal_cdf(const double z) {
    return 0.5 * std::erfc(-z / std::sqrt(2.0));
}

// Rational approximation of the standard normal quantile function
// (Acklam's algorithm), accurate to ~1.15e-9 absolute error over (0,1).
double std_normal_inv_cdf(const double p) {
    static constexpr double a[6] = {
        -3.969683028665376e+01,  2.209460984245205e+02, -2.759285104469687e+02,
         1.383577518672690e+02, -3.066479806614716e+01,  2.506628277459239e+00
    };
    static constexpr double b[5] = {
        -5.447609879822406e+01,  1.615858368580409e+02, -1.556989798598866e+02,
         6.680131188771972e+01, -1.328068155288572e+01
    };
    static constexpr double c[6] = {
        -7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00,
        -2.549732539343734e+00,  4.374664141464968e+00,  2.938163982698783e+00
    };
    static constexpr double d[4] = {
        7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00, 3.754408661907416e+00
    };
    static constexpr double p_low = 0.02425;

    err_checkf(p > 0.0 && p < 1.0, "std_normal_inv_cdf: p must be in (0,1)", std::cout);

    if (p < p_low) {
        const double q = std::sqrt(-2.0 * std::log(p));
        return (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) /
            ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0);
    }
    if (p > 1.0 - p_low) {
        const double q = std::sqrt(-2.0 * std::log(1.0 - p));
        return -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) /
            ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0);
    }
    const double q = p - 0.5;
    const double r = q * q;
    return (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q /
        (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1.0);
}

double anderson_darling_statistic(vec z) {
    const int n = static_cast<int>(z.size());
    err_checkf(n >= 2, "anderson_darling_statistic: need at least 2 values", std::cout);
    std::sort(z.begin(), z.end());

    double sum = 0.0;
    for (int i = 1; i <= n; i++) {
        // Clamp to avoid log(0) from finite-precision saturation of Phi at the tails.
        const double phi_lo = std::clamp(std_normal_cdf(z[i - 1]), 1e-300, 1.0 - 1e-16);
        const double phi_hi = std::clamp(std_normal_cdf(z[n - i]), 1e-300, 1.0 - 1e-16);
        sum += (2.0 * i - 1.0) * (std::log(phi_lo) + std::log(1.0 - phi_hi));
    }
    return -static_cast<double>(n) - sum / static_cast<double>(n);
}

ProbabilityPlotFit normal_probability_plot_fit(vec z) {
    const int n = static_cast<int>(z.size());
    err_checkf(n >= 2, "normal_probability_plot_fit: need at least 2 values", std::cout);
    std::sort(z.begin(), z.end());

    vec x(n);
    for (int i = 0; i < n; i++) {
        x[i] = std_normal_inv_cdf((static_cast<double>(i) + 0.5) / static_cast<double>(n));
    }

    double mean_x = 0.0, mean_y = 0.0;
    for (int i = 0; i < n; i++) {
        mean_x += x[i];
        mean_y += z[i];
    }
    mean_x /= n;
    mean_y /= n;

    double cov_xy = 0.0, var_x = 0.0;
    for (int i = 0; i < n; i++) {
        const double dx = x[i] - mean_x;
        cov_xy += dx * (z[i] - mean_y);
        var_x += dx * dx;
    }

    ProbabilityPlotFit fit;
    fit.slope = (var_x > 0.0) ? cov_xy / var_x : 1.0;
    fit.intercept = mean_y - fit.slope * mean_x;
    return fit;
}

double sample_skewness(const vec& z) {
    const int n = static_cast<int>(z.size());
    err_checkf(n >= 2, "sample_skewness: need at least 2 values", std::cout);
    double mean = 0.0;
    for (const double v : z) mean += v;
    mean /= n;

    double m2 = 0.0, m3 = 0.0;
    for (const double v : z) {
        const double d = v - mean;
        m2 += d * d;
        m3 += d * d * d;
    }
    m2 /= n;
    m3 /= n;
    return (m2 > 0.0) ? m3 / std::pow(m2, 1.5) : 0.0;
}

double sample_excess_kurtosis(const vec& z) {
    const int n = static_cast<int>(z.size());
    err_checkf(n >= 2, "sample_excess_kurtosis: need at least 2 values", std::cout);
    double mean = 0.0;
    for (const double v : z) mean += v;
    mean /= n;

    double m2 = 0.0, m4 = 0.0;
    for (const double v : z) {
        const double d = v - mean;
        const double d2 = d * d;
        m2 += d2;
        m4 += d2 * d2;
    }
    m2 /= n;
    m4 /= n;
    return (m2 > 0.0) ? (m4 / (m2 * m2)) - 3.0 : 0.0;
}

double jarque_bera_statistic(const vec& z) {
    const int n = static_cast<int>(z.size());
    const double s = sample_skewness(z);
    const double k = sample_excess_kurtosis(z);
    return (static_cast<double>(n) / 6.0) * (s * s + (k * k) / 4.0);
}

namespace {
    // Ranks for Spearman correlation, with average ranks for ties.
    vec rank_average(const vec& v) {
        const int n = static_cast<int>(v.size());
        ivec order(n);
        for (int i = 0; i < n; i++) order[i] = i;
        std::sort(order.begin(), order.end(), [&v](int a, int b) { return v[a] < v[b]; });

        vec ranks(n, 0.0);
        int i = 0;
        while (i < n) {
            int j = i;
            while (j + 1 < n && v[order[j + 1]] == v[order[i]]) j++;
            const double avg_rank = 0.5 * (static_cast<double>(i) + static_cast<double>(j)) + 1.0;
            for (int k = i; k <= j; k++) ranks[order[k]] = avg_rank;
            i = j + 1;
        }
        return ranks;
    }

    double spearman_correlation(const vec& a, const vec& b) {
        const vec ra = rank_average(a);
        const vec rb = rank_average(b);
        const int n = static_cast<int>(a.size());

        double mean_a = 0.0, mean_b = 0.0;
        for (int i = 0; i < n; i++) { mean_a += ra[i]; mean_b += rb[i]; }
        mean_a /= n;
        mean_b /= n;

        double cov = 0.0, var_a = 0.0, var_b = 0.0;
        for (int i = 0; i < n; i++) {
            const double da = ra[i] - mean_a;
            const double db = rb[i] - mean_b;
            cov += da * db;
            var_a += da * da;
            var_b += db * db;
        }
        return (var_a > 0.0 && var_b > 0.0) ? cov / std::sqrt(var_a * var_b) : 0.0;
    }
}

BinnedTrend binned_z_squared_trend(const vec& z, const vec& key, int n_bins, double flag_threshold) {
    const int n = static_cast<int>(z.size());
    err_checkf(z.size() == key.size(), "binned_z_squared_trend: z and key size mismatch", std::cout);
    n_bins = std::max(2, std::min(n_bins, n));

    ivec order(n);
    for (int i = 0; i < n; i++) order[i] = i;
    std::sort(order.begin(), order.end(), [&key](int a, int b) { return key[a] < key[b]; });

    vec bin_mean_z2(n_bins, 0.0);
    vec bin_center(n_bins, 0.0);
    const int base = n / n_bins;
    const int rem = n % n_bins;
    int cursor = 0;
    for (int bin = 0; bin < n_bins; bin++) {
        const int count = base + (bin < rem ? 1 : 0);
        double sum_z2 = 0.0, sum_key = 0.0;
        for (int k = 0; k < count; k++) {
            const int idx = order[cursor + k];
            sum_z2 += z[idx] * z[idx];
            sum_key += key[idx];
        }
        bin_mean_z2[bin] = (count > 0) ? sum_z2 / count : 0.0;
        bin_center[bin] = (count > 0) ? sum_key / count : 0.0;
        cursor += count;
    }

    vec bin_index(n_bins);
    for (int bin = 0; bin < n_bins; bin++) bin_index[bin] = static_cast<double>(bin);

    double mean_x = 0.0, mean_y = 0.0;
    for (int bin = 0; bin < n_bins; bin++) { mean_x += bin_index[bin]; mean_y += bin_mean_z2[bin]; }
    mean_x /= n_bins;
    mean_y /= n_bins;

    double cov_xy = 0.0, var_x = 0.0;
    for (int bin = 0; bin < n_bins; bin++) {
        const double dx = bin_index[bin] - mean_x;
        cov_xy += dx * (bin_mean_z2[bin] - mean_y);
        var_x += dx * dx;
    }

    BinnedTrend trend;
    trend.slope = (var_x > 0.0) ? cov_xy / var_x : 0.0;
    trend.spearman_r = spearman_correlation(bin_index, bin_mean_z2);
    trend.flagged = std::abs(trend.spearman_r) > flag_threshold;
    return trend;
}
