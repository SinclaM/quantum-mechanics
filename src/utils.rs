use std::fmt;

/// A constant size vector which grows inward, from both left
/// and right. 
pub struct InwardVec<T> {
    data: Vec<T>,
    left: usize,
    right: usize,
    meetup_idx: usize,
}

impl<T: Default + Clone> InwardVec<T> {
    /// Creates a new InwardVec with a given size and meetup point. The internal
    /// vector is initialized with defualt values, such that self.data.len == self.data.capacity
    /// == size.
    pub fn new(size: usize, final_left_idx: usize) -> Result<InwardVec<T>, InwardVecError> {
        if final_left_idx + 1 > size {
            Err(InwardVecError { message: String::from("Invalid meeting index.")})
        } else {
            Ok(InwardVec {
                data: vec![Default::default(); size],
                left: 0,
                right: size - 1,
                meetup_idx: final_left_idx,
            })
        }
    }

    /// Pushes a value from to the left, toward the right. This operation will fail
    /// if there is no more room in the lefthand partition.
    pub fn push_from_left(&mut self, val: T) -> Result<(), InwardVecError> {
        if self.left > self.right {
            Err(InwardVecError { message: String::from("Unable to push. Vector is full.") })
        } else if self.left > self.meetup_idx {
            Err(InwardVecError { message: String::from("Unable to push. No more rightward space.") })
        } else {
            self.data[self.left] = val; self.left += 1;
            Ok(())
        }
    }

    /// Pushes a value from to the right, toward the left. This operation will fail
    /// if there is no more room in the righthand partition.
    pub fn push_from_right(&mut self, val: T) -> Result<(), InwardVecError> {
        if self.left > self.right {
            Err(InwardVecError { message: String::from("Unable to push. Vector is full.") })
        } else if self.right < self.meetup_idx {
            Err(InwardVecError { message: String::from("Unable to push. No more leftward space.") })
        } else {
            self.data[self.right] = val;
            self.right -= 1;
            Ok(())
        }
    }

    /// Returns a reference to the internal vector 
    pub fn vec(&self) -> &Vec<T> {
        &self.data
    }

    /// Returns a mutable reference to the internal vector 
    pub fn vec_mut(&mut self) -> &mut Vec<T> {
        &mut self.data
    }

}


#[derive(Debug)]
pub struct InwardVecError {
    message: String,
}

impl fmt::Display for InwardVecError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", self.message)
    }
}

#[cfg(test)]
mod tests {
    use crate::utils::InwardVec;
    #[test]
    fn create() {
        let v: InwardVec<i32> = InwardVec::new(5, 2).unwrap();
        assert_eq!(v.vec(), &vec![0_i32; 5]);

        assert!(InwardVec::<f64>::new(0, 1).is_err());
    }

    #[test]
    fn push_until_full() {
        let mut v = InwardVec::<f64>::new(5, 2).unwrap();

        assert!(v.push_from_left(1.0).is_ok());
        assert_eq!(v.vec(), &vec![1.0, 0.0, 0.0, 0.0, 0.0]);

        assert!(v.push_from_left(2.0).is_ok());
        assert_eq!(v.vec(), &vec![1.0, 2.0, 0.0, 0.0, 0.0]);

        assert!(v.push_from_right(5.0).is_ok());
        assert_eq!(v.vec(), &vec![1.0, 2.0, 0.0, 0.0, 5.0]);

        assert!(v.push_from_left(3.0).is_ok());
        assert_eq!(v.vec(), &vec![1.0, 2.0, 3.0, 0.0, 5.0]);

        // No pushing past meetup point
        assert!(v.push_from_left(4.0).is_err());

        assert!(v.push_from_right(4.0).is_ok());
        assert_eq!(v.vec(), &vec![1.0, 2.0, 3.0, 4.0, 5.0]);

        assert!(v.push_from_right(6.0).is_err());
    }

    #[test]
    fn push_from_left() {
        let mut v = InwardVec::<i32>::new(4, 2).unwrap();
        assert!(v.push_from_left(1).is_ok());
        assert!(v.push_from_left(1).is_ok());
        assert!(v.push_from_left(1).is_ok());

        assert!(v.push_from_left(1).is_err());
    }
    
    #[test]
    fn push_from_right() {
        let mut v = InwardVec::<i32>::new(4, 2).unwrap();

        assert!(v.push_from_right(1).is_ok());
        assert!(v.push_from_right(1).is_ok());
        assert!(v.push_from_right(1).is_err());
    }
}
