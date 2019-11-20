
use lllreduce::{
		Basetype, 
        gram_schmidt_with_coeffs, 
        lll_reduce};

#[test]
fn main() {
	let original_mtx : std::vec::Vec<std::vec::Vec::<Basetype>> = vec![
		vec![0.0,3.0,4.0,7.0,8.0],
		vec![1.0,0.0,1.0,8.0,7.0],
		vec![1.0,1.0,3.0,5.0,6.0],
		vec![0.0,3.0,4.0,7.0,6.0],
		vec![0.0,3.0,4.0,8.0,9.0]
	];
    let mut mtxtuple = gram_schmidt_with_coeffs(original_mtx);
    println!("\tumtx BEFORE");
    for a in &mtxtuple.0 {
        println!("\t\t{:?}", a);
    }
    println!("\tqmtx BEFORE");
    for a in &mtxtuple.1 {
        println!("\t\t{:?}", a);
    }
    println!("\tomtx BEFORE");
    for a in &mtxtuple.2 {
        println!("\t\t{:?}", a);
    }
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
