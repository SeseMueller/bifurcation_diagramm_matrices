use ndarray_linalg::{Inverse, Norm};

const BUCKETS: usize = 1024;
const SUBBUCKETS: usize = 2048;
const ITERATIONS: usize = 100;

const MIN_X: f64 = 0.0;
const MAX_X: f64 = 1.0;

const MIN_Y: f64 = 3.0;
const MAX_Y: f64 = 4.0;

extern crate blas_src; // I have no idea why this is needed.

/// Instead of the usual, iterative approach to generating a logistic map, this program uses a
/// pregenerated map to guide in the generation. Specifically, the r-value of a specific row
/// is used to generate a lookup matrix. Specifically, the row corresponding to the r-value
/// is split into buckets. That means that each bucket represents a range of posible values
/// for the fixed r-values. The matrix then contains for each pair of buckets the fraction
/// of numbers that go from the first to the second bucket on that r-value.

/// Pregenerates a bucket vector for the given r-value or factor.
fn pregenerate(r: f64) -> ndarray::Array2<f64> {
    // For each bucket, sample the range of values that fall into it.

    let mut returnval = ndarray::Array2::zeros((BUCKETS, BUCKETS));

    for bucket in 0..BUCKETS {
        // Maps the number of that bucket to the continues range of values that fall into it.
        let min = MIN_X + (MAX_X - MIN_X) * (bucket as f64 / BUCKETS as f64);
        let max = MIN_X + (MAX_X - MIN_X) * ((bucket + 1) as f64 / BUCKETS as f64);

        //DEBUG
        // println!("Bucket {} has range {} to {}", bucket, min, max);

        // Sample the range of values that fall into the bucket.

        for s in 0..SUBBUCKETS {
            let value = min + (max - min) * (s as f64 / SUBBUCKETS as f64);
            let mapped_value = r * value * (1. - value);

            let mapped_bucket = (mapped_value - MIN_X) / (MAX_X - MIN_X) * BUCKETS as f64;

            let res_bucket = mapped_bucket.floor() as usize;

            //DEBUG
            // println!("Value {} maps to {} and to bucket {}", value, mapped_value, res_bucket);

            returnval[[res_bucket, bucket]] += 1.0 / SUBBUCKETS as f64;
        }
    }

    returnval
}

fn main() {

    println!("Adding Accelerate framework");
    println!("cargo:rustc-link-lib=framework=Accelerate");
    
    println!("Generating the logistic map");

    let test = pregenerate(3.9);

    // println!("{:?}", test);

    // debug_iterate(test);

    find_stable_point(test);

    // Paint using iteration
    bucket_iterate_and_paint();

    // Paint approximation using matrix inversion
    // bucket_inverse_and_paint();
}

// fn debug_iterate(input: ndarray::Array2<f64>) {

//     // Start off with each bucket being equally likely.
//     let mut current = ndarray::Array1::from_elem(BUCKETS, 1.0 / BUCKETS as f64);

//     // Iterate the map.
//     for _ in 0..ITERATIONS {
//         current = input.dot(&current);
//     }
//     println!("{:?}", current);

//     dbg!(input.inv());

//     // Does .inv even work on my machine?
//     let test: ndarray::Array2<f64> = ndarray::Array2::eye(3);
//     dbg!(test.inv());

// }

fn find_stable_point(buckets: ndarray::Array2<f64>) -> ndarray::Array1<f64>{
    
    // Tries to find a vector such that buckets.dot(vector) = vector. 
    // This can be theoretically multiplied by buckets.inv() to get the stable point.
    // So we are trying to solve the linear equation buckets.dot(vector) - vector = 0.
    // This can be done by using the inverse of buckets as a transformation matrix.

    // The problem is that the inverse of buckets is not always defined.
    // So we need to find a way to find a vector that is close enough to the solution.

    // As a workaround, we can use the fact that the points in our simulation *could*
    // *stay* where they are. Therefore, we could try to multiply the matrix by a value
    // such as 0.8 and add the identity matrix times 0.2. 

    // let temp = buckets.clone();
    let temp = buckets.clone() * 0.8 + ndarray::Array2::<f64>::eye(BUCKETS) * 0.2;

    // NOW:
    // for simplicity, buckets = A, vector = x, identity = I
    // Ax=x
    // Ax-x=0
    // Ax-Ix=0
    // (A-I)x=0
    // x = 0
    // .....
    // Ok
    // so
    // we could instead say that the result should only **approximate** the input, so we get
    // (A-I)x=epsilon
    // and then we can use the inverse of A-I to get
    // x = (A-I)^-1 * epsilon

    // dbg!(temp.clone() - ndarray::Array2::<f64>::eye(BUCKETS));
    // dbg!((temp.clone() - ndarray::Array2::<f64>::eye(BUCKETS)).inv().unwrap().dot(&ndarray::Array1::ones(BUCKETS)));
    // This works... sometimes? But not always. Chaos is weird.

    let inverse = (temp.clone() - ndarray::Array2::<f64>::eye(BUCKETS)).inv();

    match inverse {
        Ok(inverse) => {
            // Now we can use the inverse to get the vector (approximately).
            let vector = inverse.dot(&ndarray::Array1::ones(BUCKETS));
            // dbg!(vector.clone());
            // dbg!(buckets.dot(&vector));
            vector
        },
        Err(e) => {
            // It didn't work. As a semi-last resort, use find_almost_inverse_with_random.
            dbg!(e);
            let temp = find_almost_inverse_with_random(temp);
            let vector = temp.dot(&ndarray::Array1::ones(BUCKETS));
            // dbg!(vector.clone());
            // dbg!(buckets.dot(&vector));
            vector
        }
    }
}

fn find_almost_inverse_with_random(buckets: ndarray::Array2<f64>) -> ndarray::Array2<f64> {
    // Given an array from the above function, tries to find a matrix that is close to the inverse
    // of the input matrix. This is done by using a probabalistic method.

    // Start off with a random matrix.
    let current: ndarray::Array2<f64> = (ndarray_linalg::generate::random((BUCKETS, BUCKETS)) - 0.5) / 1000.0; //Scale from 0-1 to 0-0.001

    let inverse = (buckets.clone() + current - ndarray::Array2::<f64>::eye(BUCKETS)).inv();

    match inverse {
        Ok(inverse) => {
            // Yay! it worked!
            inverse
        }
        Err(e) => {
            // It didn't work. Try again.
            dbg!(e);
            find_almost_inverse_with_random(buckets)
        }
    }


}

fn bucket_iterate_and_paint() {
    // Each Bucket is exactly one pixel. The picture is bound by the constants defined at the top.
    // The picture is square, so the number of buckets is the same in both dimensions.

    let mut img = image::ImageBuffer::new(BUCKETS as u32, BUCKETS as u32);

    println!("Pregenerating matrices...");

    // This would be nice, but I don't have 64GB of RAM.
    // let bucket_vector_vector = (0..BUCKETS).into_iter().map(|y| {
    //     let y_value = MIN_Y + (MAX_Y - MIN_Y) * (y as f64 / BUCKETS as f64);
    //     // The y_value is the r-value for this collumn.

    //     // Pregenerate the bucket vector for this r-value.
    //     pregenerate(y_value).t().into_owned()
    // }).collect::<Vec<_>>();

    println!("Iterating and painting");

    let base_copy = ndarray::Array1::from_elem(BUCKETS, 1.0 / BUCKETS as f64);

    // Loop through BUCKET collumns, one for each pixel.
    for y in 0..BUCKETS {

        print!("Working on iteration {} of {}...\r", y, BUCKETS);

        // Pregenerate the bucket vector for this r-value.
        let bucket_vector = pregenerate(MIN_Y + (MAX_Y - MIN_Y) * (y as f64 / BUCKETS as f64));

        // Start off with each bucket being equally likely.
        let mut current = base_copy.clone();
        let mut previous = current.clone();

        // Iterate the map 
        for _ in 0..ITERATIONS {
            current = bucket_vector.dot(&current);

            // Print norm of difference between current and previous.
            let norm = (&current - &previous).norm_l1();
            // println!("{}", norm);
            if norm < 0.001 {
                // println!("Converged!");
                break;
            }
            
            previous = current.clone();
        }

        
        // Iterate the map 
        // for _ in 0..(ITERATIONS as f64).log2().ceil() as usize {
        //     bucket_vector = bucket_vector.dot(&bucket_vector);
        // }
        // current = bucket_vector.dot(&current);
        
        // Note: Iterations can be quite low. Therefore, this log2 approach is not really faster.

        // However, because we are basically trying to solve a linear equation, we could use the
        // inverse of the matrix to speed up the calculation. See bucket_inverse_and_paint.
            
        // Paint the current vector into the image.
        for x in 0..BUCKETS {
            let scaling = current[x].sqrt().sqrt();
            let value = 255.0 - (scaling * 255.0);
            img.put_pixel(y as u32, x as u32, image::Luma([value as u8]));
        }
    }

    println!();
    
    //Save the image
    let filename = format!("logistic_map_{}_{}_{}_{}_{}.png", BUCKETS, SUBBUCKETS, ITERATIONS, MIN_X, MAX_X);
    let result = img.save(filename);
    match result {
        Ok(_) => println!("Image saved successfully"),
        Err(e) => println!("Error while saving image: {}", e),
    }

}

/// Instead of iterating, works with the inverse of the matrix.
fn bucket_inverse_and_paint() {
    // Each Bucket is exactly one pixel. The picture is bound by the constants defined at the top.
    // The picture is square, so the number of buckets is the same in both dimensions.

    let mut img = image::ImageBuffer::new(BUCKETS as u32, BUCKETS as u32);

    // This would be nice, but I don't have 64GB of RAM.
    // let bucket_vector_vector = (0..BUCKETS).into_iter().map(|y| {
    //     let y_value = MIN_Y + (MAX_Y - MIN_Y) * (y as f64 / BUCKETS as f64);
    //     // The y_value is the r-value for this collumn.

    //     // Pregenerate the bucket vector for this r-value.
    //     pregenerate(y_value).t().into_owned()
    // }).collect::<Vec<_>>();

    println!("Computing and painting");

    // Loop through BUCKET collumns, one for each pixel.
    for y in 0..BUCKETS {

        print!("Working on iteration {} of {}...\r", y, BUCKETS);

        // Pregenerate the bucket vector for this r-value.
        let bucket_vector = pregenerate(MIN_Y + (MAX_Y - MIN_Y) * (y as f64 / BUCKETS as f64));

        // Start off with each bucket being equally likely.

        let stable_unnorm_vec = find_stable_point(bucket_vector);
        let stable_vec = stable_unnorm_vec.clone() / stable_unnorm_vec.norm_max();
            
        // Paint the current vector into the image.
        for x in 0..BUCKETS {
            let scaling = stable_vec[x].abs().sqrt().sqrt();
            let value = 255.0 - (scaling * 255.0);
            img.put_pixel(y as u32, x as u32, image::Luma([value as u8]));
        }
    }
    
    println!();

    //Save the image
    let filename = format!("logistic_map_{}_{}_{}_{}_{}.png", BUCKETS, SUBBUCKETS, ITERATIONS, MIN_X, MAX_X);
    let result = img.save(filename);
    match result {
        Ok(_) => println!("Image saved successfully"),
        Err(e) => println!("Error while saving image: {}", e),
    }

}