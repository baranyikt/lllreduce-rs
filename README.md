# Lenstra-Lenstra-Lovász (LLL) reduction 

Transforms a lattice's basis into a form in which the first vector of the basis is not "much" longer than the shortest (non-zero) vector of the lattice:

` ||first basis vector|| <= 2^((n-1)/2) ||shortest vector||`

where `n` is the dimension of lattice (assuming `LOVASZ_FACTOR = 4.0/3.0`).

[LLL Reduction at Wikipedia](https://en.wikipedia.org/wiki/Lenstra%E2%80%93Lenstra%E2%80%93Lov%C3%A1sz_lattice_basis_reduction_algorithm)

# Usage

```rust

use lllreduce::{
	Basetype, 
        gram_schmidt_with_coeffs, 
        lll_reduce
	};

fn main() {
	let original_mtx : std::vec::Vec<std::vec::Vec::<Basetype>> = vec![
		vec![0.0,3.0,4.0,7.0,8.0],
		vec![1.0,0.0,1.0,8.0,7.0],
		vec![1.0,1.0,3.0,5.0,6.0],
		vec![0.0,3.0,4.0,7.0,6.0],
		vec![0.0,3.0,4.0,8.0,9.0]
	];
    let mut mtxtuple = lllreduce::gram_schmidt_with_coeffs(original_mtx);
	lll_reduce(&mut mtxtuple);
    println!("\tLLL-reduced basis");
    for a in &mtxtuple.2 {
        println!("\t\t{:?}", a);
    }    
    println!("\tumtx");
    for a in &mtxtuple.0 {
        println!("\t\t{:?}", a);
    }
    println!("\tqmtx");
    for a in &mtxtuple.1 {
        println!("\t\t{:?}", a);
    }
}
```

# Note
This is the very first, very drafty version, anything can change in it in the future.
