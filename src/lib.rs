use wasm_bindgen::prelude::*;

#[wasm_bindgen]
pub fn calculate_vincenty(lat1: f64, lon1: f64, lat2: f64, lon2: f64) -> f64 {
    const A: f64 = 6378137.0; // Semi-major axis (meters)
    const F: f64 = 1.0 / 298.257223563; // Flattening
    const B: f64 = A * (1.0 - F); // Semi-minor axis

    let phi1 = lat1.to_radians();
    let phi2 = lat2.to_radians();
    let lambda1 = lon1.to_radians();
    let lambda2 = lon2.to_radians();

    let u1 = ((1.0 - F) * phi1.tan()).atan();
    let u2 = ((1.0 - F) * phi2.tan()).atan();

    let mut lambda = lambda2 - lambda1;
    let mut sin_sigma;
    let mut cos_sigma;
    let mut sigma;
    let mut sin_alpha;
    let mut cos2_alpha;
    let mut cos2_sigma_m;
    let mut lambda_prev;
    const MAX_ITERATIONS: u32 = 1000;
    const CONVERGENCE_THRESHOLD: f64 = 1e-12;

    let mut iteration_count = 0;

    loop {
        let sin_lambda = lambda.sin();
        let cos_lambda = lambda.cos();
        sin_sigma = ((u2.cos() * sin_lambda).powi(2)
            + (u1.cos() * u2.sin() - u1.sin() * u2.cos() * cos_lambda).powi(2))
            .sqrt();
        cos_sigma = u1.sin() * u2.sin() + u1.cos() * u2.cos() * cos_lambda;
        sigma = sin_sigma.atan2(cos_sigma);
        sin_alpha = u1.cos() * u2.cos() * sin_lambda / sin_sigma;
        cos2_alpha = 1.0 - sin_alpha.powi(2);
        cos2_sigma_m = if cos2_alpha != 0.0 {
            cos_sigma - 2.0 * u1.sin() * u2.sin() / cos2_alpha
        } else {
            0.0
        };

        let c = F / 16.0 * cos2_alpha * (4.0 + F * (4.0 - 3.0 * cos2_alpha));
        lambda_prev = lambda;
        lambda = (lambda2 - lambda1)
            + (1.0 - c)
                * F
                * sin_alpha
                * (sigma
                    + c * sin_sigma * (cos2_sigma_m + c * cos_sigma * (-1.0 + 2.0 * cos2_sigma_m.powi(2))));

        iteration_count += 1;
        if (lambda - lambda_prev).abs() < CONVERGENCE_THRESHOLD || iteration_count > MAX_ITERATIONS {
            break;
        }
    }

    let u2 = cos2_alpha * (A.powi(2) - B.powi(2)) / B.powi(2);
    let a = 1.0 + u2 / 16384.0 * (4096.0 + u2 * (-768.0 + u2 * (320.0 - 175.0 * u2)));
    let b = u2 / 1024.0 * (256.0 + u2 * (-128.0 + u2 * (74.0 - 47.0 * u2)));
    let delta_sigma = b
        * sin_sigma
        * (cos2_sigma_m + b / 4.0 * (cos_sigma * (-1.0 + 2.0 * cos2_sigma_m.powi(2))
            - b / 6.0 * cos2_sigma_m * (-3.0 + 4.0 * sin_sigma.powi(2)) * (-3.0 + 4.0 * cos2_sigma_m.powi(2))));
    B * a * (sigma - delta_sigma)
}
