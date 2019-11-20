//! A libaray for running Lenstra-Lenstra-Lovasz (LLL) reduction on grid bases

use std::mem;

const LOVASZ_FACTOR : Basetype = 4.0/3.0;
const LOG_STEPS : bool = false;
pub type Basetype = f32;
pub type MtxOfBaseType = std::vec::Vec<std::vec::Vec<Basetype>>;
pub type MatrixTriplet = (MtxOfBaseType, MtxOfBaseType, MtxOfBaseType);

pub trait DotProduct<BaseField> {
	type Output;
	fn dotprod(&self, other: &Self) -> BaseField;
}

impl DotProduct<Basetype> for std::vec::Vec<Basetype> {
	type Output = Basetype;
	fn dotprod(&self, other: &Self) -> Basetype {
		let mut retval: Basetype = 0.0;
		for it in self.iter().zip(other.iter()) {
			let (a,b) = it;
			retval += a * b;
		}
		retval
	}
}

fn multiply_by_scalar(v: &std::vec::Vec<Basetype>, s: Basetype) -> std::vec::Vec<Basetype>{
	let mut retval : std::vec::Vec<Basetype> = v.clone();
	for elem in &mut retval {
		*elem *= s;
	}
	retval
}


fn add_vector(v: &mut std::vec::Vec<Basetype>, a: &std::vec::Vec<Basetype>) {
	for it in v.iter_mut().zip(a.iter()) {
		let (vi, ai) = it;
		*vi += ai;
	}
}

fn sub_vector(v: &mut std::vec::Vec<Basetype>, a: &std::vec::Vec<Basetype>) {
	for it in v.iter_mut().zip(a.iter()) {
		let (vi, ai) = it;
		*vi -= ai;
	}
}

pub fn gram_schmidt_with_coeffs(mtx: MtxOfBaseType) -> MatrixTriplet {		
	let mut umtx: MtxOfBaseType = MtxOfBaseType::with_capacity(mtx.len());
	let mut qmtx: MtxOfBaseType = MtxOfBaseType::with_capacity(mtx.len());
	for vec_ind in 0..mtx.len() {
		let mut u = mtx[vec_ind].clone(); // std::vec::Vec::<Basetype>::with_capacity(vec.len());
		let mut q = std::vec::Vec::<Basetype>::with_capacity(u.len());
		for orto_ind in 0..vec_ind {
			let coeff = mtx[vec_ind].dotprod(&umtx[orto_ind])/umtx[orto_ind].dotprod(&umtx[orto_ind]);
			sub_vector(&mut u, &multiply_by_scalar(&umtx[orto_ind], coeff));
			q.push(coeff);
		}
		q.push(1.0);
		for _ in vec_ind+1..u.len() {
			q.push(0.0);
		}
		qmtx.push(q);
		umtx.push(u);
	}
	(umtx, qmtx, mtx)
}

fn gram_schmidt_with_coeffs_update_at_rowpair(mtxtuple: &mut MatrixTriplet, i: usize) {
	let (umtx, qmtx, origmtx) = mtxtuple;
	for vec_ind in i..i+2 {
		let mut u = origmtx[vec_ind].clone();
		let q = &mut qmtx[vec_ind];
		for orto_ind in 0..vec_ind {
			let coeff = origmtx[vec_ind].dotprod(&umtx[orto_ind])/umtx[orto_ind].dotprod(&umtx[orto_ind]);
			sub_vector(&mut u, &multiply_by_scalar(&umtx[orto_ind], coeff));
			q[orto_ind] = coeff;
		}
		q[vec_ind] = 1.0;						//TODO maybe superfluous
		for orto_ind in vec_ind+1..u.len() {
			q[orto_ind] = 0.0;
		}
		mem::swap(&mut umtx[vec_ind], &mut u);
	}
}


/// function primitive 1 for MatrixTriplet: adding lambda * v[j] to v[i] 		
/// [prereq: 0 <= j < i < mtx.len()]
fn mtxtuple_primitive1(mtxtuple: &mut MatrixTriplet, i: usize, j: usize, lambda: Basetype) {		
	let (_umtx, qmtx, omtx) = mtxtuple;
	let (ohead, otail) = omtx.split_at_mut(j + 1);
	add_vector(&mut ohead[j], &multiply_by_scalar(&otail[i - j - 1], lambda));
	for t in 0..j+1 {
		if LOG_STEPS {
			println!("qmtx[{}][{}] += lambda * qmtx[{}][{}]", i, t, j, t);
		}
		qmtx[i][t] += lambda * qmtx[j][t];
	}
}

fn gauss_reducer(mtxtuple: &mut MatrixTriplet) -> bool {		
	let mut did_anything = false;
	let mtxlen = mtxtuple.0.len();
	for j in (0..mtxlen).rev() {
		for i in j+1..mtxlen {
			let is_bad_i_j = mtxtuple.1[i][j] > 0.5;
			let rounded = mtxtuple.1[i][j].round();
			if LOG_STEPS {
				println!("{},{}\t{}\t{:?}\t{}", i,j, mtxtuple.1[i][j], is_bad_i_j, rounded);
			}
			if is_bad_i_j {
				mtxtuple_primitive1(mtxtuple, i, j, -rounded);
				did_anything = true;
			}
		}
	}
	did_anything
}

fn mtxtuple_primitive2(mtxtuple: &mut MatrixTriplet, i: usize) {		
	let (_umtx, _qmtx, omtx) = mtxtuple;
	let (head, tail) = omtx.split_at_mut(i + 1);
	mem::swap(&mut head[i], &mut tail[0]);
	gram_schmidt_with_coeffs_update_at_rowpair(mtxtuple, i);
}

fn calculate_orthogonal_component_len_sqr(mtxtuple: &MatrixTriplet, vector_index: usize, subspace_index: usize) -> Basetype {		
	let (umtx, qmtx, omtx) = mtxtuple;
	let mut orthogonal_vector = omtx[vector_index].clone();
	for j in 0..subspace_index {
		sub_vector(&mut orthogonal_vector, &multiply_by_scalar(&umtx[j], qmtx[vector_index][j]));			//TODO check if i,j or j,i
	}
	orthogonal_vector.dotprod(&orthogonal_vector)
}

fn is_blocking_lovasz_optimality(mtxtuple: &MatrixTriplet, i: usize) -> bool {		
	calculate_orthogonal_component_len_sqr(mtxtuple, i, i) > LOVASZ_FACTOR * calculate_orthogonal_component_len_sqr(mtxtuple, i+1, i)
}

/// LLL reduction step: swaps any (i,i+1) base vector pair blocking LLL-optimality (i.e. those having b' components orthogonal
/// to space V_(i-1) that satisfy the inequality ||b'_i|| > LOVASZ_FACTOR ||b'_(i+1)|| )
fn lll_reduce_step(mtxtuple: &mut MatrixTriplet) -> bool {		
	let mut did_anything = false;
	let mtxlen = mtxtuple.0.len();
	for i in 0..mtxlen-1 {
		if is_blocking_lovasz_optimality(mtxtuple, i) {
			if LOG_STEPS {
				println!("pair {}/{} is blocking lovasz optimality, swapping", i, i+1);
			}
			mtxtuple_primitive2(mtxtuple, i);
			did_anything = true;
			break;
		}
	}
	did_anything
}

/// Applies LLL-reduction to mtxtuple by iterating through gauss_reducer step and lll_reduce_step until either make any change
pub fn lll_reduce(mtxtuple: &mut MatrixTriplet) {		
	let mut gauss_step_changed_base = true;
	let mut lll_step_changed_base = true;
	let mut counter = 0;
	while gauss_step_changed_base || lll_step_changed_base {
		gauss_step_changed_base = gauss_reducer(mtxtuple);
		{
			if LOG_STEPS {
				println!("\t{}/gauss", counter);
				for a in &mtxtuple.1 {
					println!("\t\t{:?}", a);
				}
			}
		}
		lll_step_changed_base = lll_reduce_step(mtxtuple);
		{
			if LOG_STEPS {
				println!("\t{}/lll", counter);
				for a in &mtxtuple.1 {
					println!("\t\t{:?}", a);
				}
			}
		}
		counter += 1;
	}
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_gram_schmidt() {
        let original_mtx : std::vec::Vec<std::vec::Vec::<Basetype>> = vec![
            vec![0.0,3.0,4.0,7.0,8.0],
            vec![1.0,0.0,1.0,8.0,7.0],
            vec![1.0,1.0,3.0,5.0,6.0],
            vec![0.0,3.0,4.0,7.0,6.0],
            vec![0.0,3.0,4.0,8.0,9.0]
        ];
	    let (orto_u, _orto_q, _omtx) = gram_schmidt_with_coeffs(original_mtx); 
		let test_orto_u : std::vec::Vec<std::vec::Vec::<Basetype>> = vec![
			vec![0.0,			3.0,			4.0,			7.0,			8.0],
			vec![1.0,			-58.0/23.0,		-163.0/69.0,	146.0/69.0,		19.0/69.0],
			vec![957.0/1207.0,	-734.0/1207.0,	783.0/1207.0,	-494.0/1207.0,	316.0/1207.0],
			vec![94.0/333.0,	76.0/999.0,		22.0/37.0,		748.0/999.0,	-980.0/999.0],
			vec![-8.0/245.0,	-13.0/490.0,	4.0/245.0,		1.0/490.0,		0.0]
		];
        const TOLERANCE : Basetype = 10e-6;
		for o_it in orto_u.iter().zip(test_orto_u.iter()) {
			let (row,testrow) = o_it;
            for i_it in row.iter().zip(testrow.iter()) {
                let (element, test_element) = i_it;
                if (element - test_element).abs() > TOLERANCE {
                    panic!("too far");
                }
            }
        }
    }
}
