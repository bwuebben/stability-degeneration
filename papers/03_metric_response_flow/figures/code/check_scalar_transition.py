#!/usr/bin/env python3
"""Finite scalar identities and controls; not a geometric PDE verification."""

from __future__ import annotations

from decimal import Decimal, localcontext
import json
from pathlib import Path


D = Decimal
PRECISION = 240
ONE = D(1)
TWO = D(2)


def one_minus_exp_neg(x: Decimal) -> Decimal:
    assert x >= 0
    if not x:
        return D(0)
    if x >= D(1000):
        return ONE
    if x < D("0.5"):
        term = total = x
        n = 2
        while True:
            term *= -x / D(n)
            nxt = total + term
            if nxt == total:
                return total
            total = nxt
            n += 1
    return ONE - (-x).exp()


def exp_neg(x: Decimal) -> Decimal:
    assert x >= 0
    return D(0) if x >= D(1000) else (-x).exp()


def neg_log_one_minus(u: Decimal) -> Decimal:
    assert 0 <= u < 1
    if not u:
        return D(0)
    if u < D("0.5"):
        power = total = u
        n = 2
        while True:
            power *= u
            nxt = total + power / D(n)
            if nxt == total:
                return total
            total = nxt
            n += 1
    return -(ONE - u).ln()


def log_one_plus(u: Decimal) -> Decimal:
    assert u >= 0
    if not u:
        return D(0)
    if u < D("0.5"):
        power = total = u
        n = 2
        while True:
            power *= -u
            nxt = total + power / D(n)
            if nxt == total:
                return total
            total = nxt
            n += 1
    return (ONE + u).ln()


def negative_logistic_log(d: Decimal, z0: Decimal, rate: Decimal, tau: Decimal, restart: Decimal = ONE) -> Decimal:
    """Log of the negative-gap solution, valid beyond exponential underflow."""
    assert d > 0 and z0 > 0 and rate > 0 and tau >= restart
    exponent = TWO * rate * d * (tau - restart) + log_one_plus(d / z0)
    return d.ln() - exponent - one_minus_exp_neg(exponent).ln()


def negative_profile_log(d: Decimal, tau: Decimal) -> Decimal:
    exponent = TWO * d * tau
    return d.ln() - exponent - one_minus_exp_neg(exponent).ln()


def negative_unit_log_ratio(d: Decimal, z0: Decimal, tau: Decimal, restart: Decimal = ONE) -> Decimal:
    """Log(Y_1/P_-d), retaining the small offset at enormous times."""
    exponent = TWO * d * tau
    displacement = log_one_plus(d / z0) - TWO * d * restart
    shifted = exponent + displacement
    assert shifted > 0
    return -displacement + one_minus_exp_neg(exponent).ln() - one_minus_exp_neg(shifted).ln()


def profile(delta: Decimal, tau: Decimal) -> Decimal:
    assert delta >= 0 and tau > 0
    return ONE / (TWO * tau) if not delta else delta / one_minus_exp_neg(TWO * delta * tau)


def logistic(delta: Decimal, z0: Decimal, rate: Decimal, tau: Decimal, restart: Decimal = ONE) -> Decimal:
    assert delta >= 0 and z0 >= delta and z0 > 0 and rate > 0 and tau >= restart
    elapsed = tau - restart
    if not delta:
        return z0 / (ONE + TWO * rate * z0 * elapsed)
    x = TWO * rate * delta * elapsed
    denominator = one_minus_exp_neg(x) + (delta / z0) * exp_neg(x)
    return delta / denominator


def offset(delta: Decimal, z0: Decimal, rate: Decimal, restart: Decimal = ONE) -> Decimal:
    assert delta < z0
    if not delta:
        return ONE / (TWO * rate * z0) - restart
    return -restart + neg_log_one_minus(delta / z0) / (TWO * rate * delta)


def logistic_from_offset(delta: Decimal, z0: Decimal, rate: Decimal, tau: Decimal) -> Decimal:
    shifted = tau + offset(delta, z0, rate)
    assert shifted > 0
    return ONE / (TWO * rate * shifted) if not delta else delta / one_minus_exp_neg(TWO * rate * delta * shifted)


def primitive(z: Decimal, delta: Decimal, alpha: Decimal) -> Decimal:
    """Antiderivative of 1/[z(z-delta)(1+alpha*z)]."""
    assert z > max(delta, D(0)) and ONE + alpha * z > 0
    if not delta:
        return -alpha * z.ln() - ONE / z + alpha * (ONE + alpha * z).ln()
    denominator = ONE + alpha * delta
    return -z.ln() / delta + (z - delta).ln() / (delta * denominator) + alpha * (ONE + alpha * z).ln() / denominator


def run() -> dict:
    with localcontext() as ctx:
        ctx.prec = PRECISION
        tolerance = D("1e-180")
        parameter_tolerance = D("1e-130")
        values = [ONE, ONE - D("0.001"), ONE + D("0.001"), ONE - ONE / 64, ONE + ONE / 64, ONE - ONE / 16, ONE + ONE / 16]
        deltas = [D(0), D("1e-80"), D("1e-40"), D("1e-10"), D("0.0001"), D("0.01"), ONE / 32]
        ordinary_times = [ONE, D("1.0001"), D("1.1"), D(2), D(10), D("1e3"), D("1e20"), D("1e80"), D("1e100")]
        cases = forms = offset_bounds = 0
        max_error_over_eta = D(0)
        max_offset_over_eta = D(0)
        for v in values:
            for rate in values:
                z0 = v / TWO
                for delta in deltas:
                    eta = max(abs(v - ONE), abs(rate - ONE), TWO * delta)
                    assert eta <= ONE / 16 and delta < z0
                    theta = offset(delta, z0, rate)
                    assert abs(theta) <= 4 * eta + tolerance
                    if eta:
                        max_offset_over_eta = max(max_offset_over_eta, abs(theta) / eta)
                    offset_bounds += 1
                    times = set(ordinary_times)
                    if delta:
                        times.update(t / delta for t in [D("0.0001"), D("0.1"), ONE, D(10), D("1e3")] if t / delta >= 1)
                    for tau in sorted(times):
                        y = logistic(delta, z0, rate, tau)
                        p = profile(delta, tau)
                        y_offset = logistic_from_offset(delta, z0, rate, tau)
                        assert abs(y_offset / y - ONE) < tolerance, (v, rate, delta, tau)
                        forms += 1
                        err = abs(y / p - ONE)
                        assert err <= 10 * eta + tolerance, (v, rate, delta, tau, err, eta)
                        if eta:
                            max_error_over_eta = max(max_error_over_eta, err / eta)
                        cases += 1

        variable_cases = 0
        for alpha in [D(-1), ONE, TWO]:
            for delta in [D(0), D("1e-80"), D("1e-6"), D("0.001")]:
                for z0 in [D("0.01"), D("0.02")]:
                    epsilon = abs(alpha) * z0
                    for fraction in [D("0.9"), D("0.5"), D("0.1"), D("0.01"), D("1e-20")]:
                        z = delta + fraction * (z0 - delta)
                        elapsed = (primitive(z0, delta, alpha) - primitive(z, delta, alpha)) / TWO
                        assert elapsed > 0
                        tau = ONE + elapsed
                        lower = logistic(delta, z0, ONE + epsilon, tau)
                        upper = logistic(delta, z0, ONE - epsilon, tau)
                        assert lower <= z * (ONE + parameter_tolerance)
                        assert upper >= z * (ONE - parameter_tolerance)
                        variable_cases += 1

        negative_profile_cases = 0
        negative_parameter_cases = 0
        max_negative_log_error_over_eta = D(0)
        for v in values:
            z0 = v / TWO
            for d in deltas[1:]:
                eta = max(abs(v - ONE), TWO * d)
                theta = -ONE + log_one_plus(d / z0) / (TWO * d)
                assert abs(theta) <= 3 * eta + tolerance
                negative_parameter_cases += 1
                times = set(ordinary_times)
                times.update(t / d for t in [D("0.0001"), D("0.1"), ONE, D(10), D("1e3")] if t / d >= 1)
                for tau in sorted(times):
                    log_ratio = negative_unit_log_ratio(d, z0, tau)
                    assert abs(log_ratio) <= 4 * eta + tolerance
                    if eta:
                        max_negative_log_error_over_eta = max(max_negative_log_error_over_eta, abs(log_ratio) / eta)
                    # The two exponent representations agree independently of
                    # whether their exponentials can be represented directly.
                    first = TWO * d * (tau - ONE) + log_one_plus(d / z0)
                    second = TWO * d * (tau + theta)
                    assert abs(first - second) <= tolerance * max(ONE, abs(first))
                    negative_profile_cases += 1

        negative_variable_cases = 0
        max_negative_variable_log_error_over_epsilon = D(0)
        for alpha in [D(-1), ONE, TWO]:
            for d in [D("1e-80"), D("1e-6"), D("0.0001")]:
                for z0 in [D("0.01"), D("0.02")]:
                    restart = ONE / (TWO * z0)
                    epsilon = TWO * abs(alpha) * z0
                    eta = TWO * d * restart
                    assert epsilon <= D("0.5") and eta <= ONE / 16
                    for fraction in [D("0.9"), D("0.5"), D("0.1"), D("0.01"), D("1e-20"), D("1e-80"), D("1e-160")]:
                        z = fraction * z0
                        elapsed = (primitive(z0, -d, alpha) - primitive(z, -d, alpha)) / TWO
                        assert elapsed > 0
                        tau = restart + elapsed
                        exponent = TWO * d * elapsed + log_one_plus(d / z0)
                        Q = log_one_plus(d / z)
                        E = Q - exponent
                        assert abs(E) <= epsilon * min(exponent, ONE) + parameter_tolerance
                        log_model = negative_logistic_log(d, z0, ONE, tau, restart)
                        log_error = abs(z.ln() - log_model)
                        assert log_error <= 3 * epsilon + parameter_tolerance
                        max_negative_variable_log_error_over_epsilon = max(max_negative_variable_log_error_over_epsilon, log_error / epsilon)
                        log_profile_error = abs(z.ln() - negative_profile_log(d, tau))
                        assert log_profile_error <= 3 * epsilon + 4 * eta + parameter_tolerance
                        negative_variable_cases += 1

        coordinate_cases = 0
        for c in [D("0.5"), ONE, TWO]:
            for d in [D(-1), ONE]:
                for w0 in [D("0.01"), D("0.03")]:
                    epsilon = abs(d) * w0 / (c - abs(d) * w0)
                    for fraction in [D("0.001"), D("0.1"), D("0.5"), ONE]:
                        w = fraction * w0
                        z = c * w + d * w * w
                        speed = w * (c + 2 * d * w) / z
                        assert abs(speed - ONE) <= epsilon + parameter_tolerance
                        coordinate_cases += 1

        delta = D("1e-40")
        tau = ONE
        wrong_factor_profile = delta / one_minus_exp_neg(delta * tau)
        correct = profile(delta, tau)
        wrong_factor_detected = abs(correct / wrong_factor_profile - ONE) > D("0.49")
        assert wrong_factor_detected

        z0 = ONE
        restart_failure = abs(logistic(D(0), z0, ONE, ONE) / profile(D(0), ONE) - ONE)
        assert restart_failure == ONE

        delta, z0 = D("1e-6"), D("0.5")
        wrong_theta = ONE + neg_log_one_minus(delta / z0) / (TWO * delta)
        wrong_y = delta / one_minus_exp_neg(TWO * delta * (ONE + wrong_theta))
        wrong_offset_detected = abs(wrong_y / z0 - ONE) > D("0.6")
        assert wrong_offset_detected

        negative_rate_log_error = abs(
            negative_logistic_log(D("0.01"), D("0.5"), D("1.01"), D("1e10"))
            - negative_profile_log(D("0.01"), D("1e10"))
        )
        assert negative_rate_log_error > D("1e5")

        return {
            "scope": "Finite scalar identities, uniform-bound samples and failure controls; no geometric PDE verification",
            "decimal_precision": PRECISION,
            "restart_parameter_cases": offset_bounds,
            "direct_offset_identity_cases": forms,
            "uniform_profile_cases": cases,
            "smallest_positive_delta_times_restart": "1e-80",
            "largest_time_over_restart": "1e100",
            "maximum_observed_profile_error_over_eta": str(max_error_over_eta),
            "proved_profile_error_constant": 10,
            "maximum_observed_offset_over_eta_times_restart": str(max_offset_over_eta),
            "proved_offset_constant": 4,
            "variable_rate_separation_cases": variable_cases,
            "calibrated_coordinate_cases": coordinate_cases,
            "negative_restart_parameter_cases": negative_parameter_cases,
            "negative_uniform_profile_cases": negative_profile_cases,
            "negative_variable_rate_separation_cases": negative_variable_cases,
            "maximum_negative_unit_rate_log_error_over_eta": str(max_negative_log_error_over_eta),
            "proved_negative_unit_rate_log_constant": 4,
            "maximum_negative_variable_rate_log_error_over_epsilon": str(max_negative_variable_log_error_over_epsilon),
            "proved_negative_variable_rate_log_constant": 3,
            "wrong_factor_two_detected": wrong_factor_detected,
            "uncalibrated_restart_error_at_left_endpoint": str(restart_failure),
            "wrong_restart_offset_sign_detected": wrong_offset_detected,
            "negative_constant_rate_sandwich_failure_detected": True,
            "negative_constant_rate_log_error_at_time_1e10": str(negative_rate_log_error),
        }


if __name__ == "__main__":
    result = run()
    output = (Path(__file__).resolve().parents[1] / "data/scalar_transition_checks.json")
    output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
