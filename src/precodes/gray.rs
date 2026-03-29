// Copyright (c) 2025 Shenghao Yang. All rights reserved.
// Licensed under the MIT License. See LICENSE-MIT for details.

/// Generates balanced Gray code sequences of constant Hamming weight.
///
/// Enumerates all binary vectors of a given `length` with exactly `weight` ones,
/// changing exactly two bit positions between consecutive columns. Used to
/// construct HDPC binary matrices efficiently.
pub struct Gray {
    length: usize,
    weight: usize,
    current_column: Vec<bool>,
    forward_index: Vec<usize>,
    backward_index: Vec<usize>,
    forward_weight: usize,
    backward_weight: usize,
    column_num: usize,
}

impl Gray {
    /// Creates a new Gray code iterator starting at the lexicographically first column.
    pub fn new(length: usize, weight: usize) -> Self {
        let mut g: Vec<bool> = vec![true; weight];
        g.extend(vec![false; length-weight]);
        let mut forward_index: Vec<usize> = (1..=length).collect();
        forward_index[0] = weight;
        let backward_index: Vec<usize> = vec![length; length];

        Self {
            length,
            weight,
            current_column: g,
            forward_index,
            backward_index,
            forward_weight: weight,
            backward_weight: weight,
            column_num: 0,
        }
    }

    /// Resets the iterator and advances to column `col_num` by stepping forward.
    pub fn with_column(&mut self, col_num: usize) {
        self.current_column = vec![true; self.weight];
        self.current_column.extend(vec![false; self.length-self.weight]);
        self.forward_index = (1..=self.length).collect();
        self.forward_index[0] = self.weight;
        self.backward_index = vec![self.length; self.length];
        self.forward_weight = self.weight;
        self.backward_weight = self.weight;
        self.column_num = 0;
        for _ in 0..col_num {
            self.next_column();
        }
    }

    /// Moves to the previous column and returns the indices of the two flipped bits.
    pub fn previous_delta(&mut self) -> Vec<usize> {
        //if self.column_num == 0 {
        //    return (0..self.weight).collect();
        //}
        let current_column = self.current_column.clone();
        self.previous_column();
        self.current_column
            .iter()
            .zip(current_column.iter())
            .map(|(x, y)| *x ^ *y)
            .enumerate()
            .filter(|(_, x)| *x)
            .map(|(i, _)| i)
            .collect()
    }

    /// Returns the current column as a `0`/`1` byte vector.
    #[allow(dead_code)] // used by tests and debugging; not used by `R10HDPC` production path
    pub fn current_column(&self) -> Vec<u8> {
        self.current_column.iter().map(|x| if *x {1} else {0}).collect()
    }

    /// Returns the indices of the `1`-positions in the current column.
    pub fn current_column_positions(&self) -> Vec<usize> {
        self.current_column.iter().enumerate().filter(|(_, x)| **x).map(|(i, _)| i).collect()
    }

    /// Returns `(forward_index[0], forward_weight)` for the current forward traversal state.
    #[allow(dead_code)]
    pub fn current_forward_state(&self) -> (usize, usize) {
        (self.forward_index[0], self.forward_weight)
    }

    /// Returns `(backward_index[0], backward_weight)` for the current backward traversal state.
    #[allow(dead_code)]
    pub fn current_backward_state(&self) -> (usize, usize) {
        (self.backward_index[0], self.backward_weight)
    }

    /// Advances to the next column in the Gray code sequence.
    pub fn next_column(&mut self) {
        let idx = self.forward_index[0]; // the index of the top element

        if idx == self.length { // last column
            self.forward_index[0] = self.weight;
            self.current_column[self.weight-1] = true;
            self.current_column[self.length-1] = false;
            self.forward_weight = self.weight;

            self.backward_index = vec![self.length; self.length];
            self.column_num = 0;
            return;
        }

        self.forward_index[0] = self.forward_index[idx];
        self.forward_index[idx] = idx + 1;

        // update the current column
        if self.current_column[idx] {
            if self.forward_weight > 0 {
                self.current_column[self.forward_weight-1] = !self.current_column[self.forward_weight-1];
            } else {
                self.current_column[idx-1] = !self.current_column[idx-1];
            }
            self.forward_weight += 1;
        } else { // current_column[idx] == false 
            if self.forward_weight > 1 {
                self.current_column[self.forward_weight-2] = !self.current_column[self.forward_weight-2];
            } else {
                self.current_column[idx-1] = !self.current_column[idx-1];
            }
            self.forward_weight -= 1;
        }
        self.current_column[idx] = !self.current_column[idx];

        // update the backward weight
        self.backward_weight = self.forward_weight;
        /* 
        if idx < self.backward_index[0] {
            self.backward_index[idx] = self.backward_index[0];
        } else if idx < self.backward_index[self.backward_index[0]] {
            self.backward_index[idx] = self.backward_index[self.backward_index[0]];
        }*/
        let mut i = 0;
        while i < self.length {
            if idx < self.backward_index[i] {
                self.backward_index[idx] = self.backward_index[i];
                break;
            } else {
                i = self.backward_index[i];
            }
        }
        self.backward_index[0] = idx;

        // forward index and weight
        if self.forward_weight == idx || self.forward_weight == 0 {
            self.forward_weight += 1;
        } else {
            self.forward_weight -= self.current_column[idx-1] as usize;
            self.forward_index[idx-1] = self.forward_index[0];
            if self.forward_weight == 0 {
                self.forward_index[0] = idx-1;
            } else {
                self.forward_index[0] = self.forward_weight;
            }
        }
        self.column_num += 1;
    }

    /// Steps back to the previous column in the Gray code sequence.
    pub fn previous_column(&mut self) {
        let idx = self.backward_index[0]; // the index of the top element
        
        if idx == self.length { // last column
            self.backward_index = (1..=self.length).collect();
            if self.weight == 1 {
                self.backward_index[0] = self.length-1  ;
            } else {
                self.backward_index[0] = self.weight-1;
            }
            self.current_column[self.weight-1] = false;
            self.current_column[self.length-1] = true;
            self.backward_weight = self.weight-1;

            self.forward_index[0] = self.length;
            self.column_num = 0;
            return;
        }

        self.backward_index[0] = self.backward_index[idx];
        self.backward_index[idx] = idx + 1;

        // update the current column
        if self.current_column[idx] {
            if self.backward_weight > 0 {
                self.current_column[self.backward_weight-1] = !self.current_column[self.backward_weight-1];
            } else {
                self.current_column[idx-1] = !self.current_column[idx-1];
            }
            self.backward_weight += 1;
        } else { // current_column[idx] == false 
            if self.backward_weight > 1 {
                self.current_column[self.backward_weight-2] = !self.current_column[self.backward_weight-2];
            } else {
                self.current_column[idx-1] = !self.current_column[idx-1];
            }
            self.backward_weight -= 1;
        }
        self.current_column[idx] = !self.current_column[idx];

        // update the forward weight
        self.forward_weight = self.backward_weight;
        if idx < self.forward_index[0] {
            self.forward_index[idx] = self.forward_index[0];
        } else if idx < self.forward_index[idx-1] {
            self.forward_index[idx] = self.forward_index[idx-1];
        } 
        self.forward_index[0] = idx;

        // backward index and weight
        if self.backward_weight == idx || self.backward_weight == 0 {
            self.backward_weight += 1;
        } else {
            self.backward_weight -= self.current_column[idx-1] as usize;
            self.backward_index[idx-1] = self.backward_index[0];
            if self.backward_weight == 0 {
                self.backward_index[0] = idx-1;
            } else {
                self.backward_index[0] = self.backward_weight;
            }
        }
        self.column_num -= 1;
    }
}


#[cfg(test)]
mod tests {
    use super::*;

    fn binomial(n: usize) -> usize {
        if n == 1 {
            2
        } else {
            binomial(n - 1) * (4 * n - 2) / n
        }
    }

    #[test]
    fn test_gray_code_with_column() {
        let length = 12;
        let weight = length / 2;
        let mut gray = Gray::new(length, weight);
        let col_num = binomial(weight);
        println!("col_num: {}", col_num);
        gray.with_column(col_num-1);
        assert_eq!(gray.column_num, col_num-1);
        for _ in 0..col_num - 1 {
            //println!("column number: {:?}", gray.column_num);
            assert_eq!(gray.current_column().iter().filter(|&&x| x == 1).count(), weight);
            assert_eq!(gray.previous_delta().len(), 2);
        }
        assert_eq!(gray.current_column_positions().len(), weight);
    }

    #[test]
    fn test_gray_code_backward() {
        let mut gray = Gray::new(5, 2);

        for i in 0..9 {
            println!("i: {}", i);
            gray.next_column();
            println!("current_column: {:?}", gray.current_column());
            println!("current state: forward {:?}, backward {:?}", gray.current_forward_state(), gray.current_backward_state());
            println!("forward_index: {:?}", gray.forward_index);
            println!("backward_index: {:?}", gray.backward_index);
        }

        for i in 0..9 {
            println!("i: {}", i);
            gray.previous_column();
            println!("current_column: {:?}", gray.current_column());
            println!("current state: forward {:?}, backward {:?}", gray.current_forward_state(), gray.current_backward_state());
            println!("forward_index: {:?}", gray.forward_index);
            println!("backward_index: {:?}", gray.backward_index);
        }
    }

    #[test]
    fn test_gray_code_5_2_forward() {
        let mut gray = Gray::new(5, 2);
        let state = [
            (2, 2),
            (1, 0),
            (3, 2),
            (2, 0),
            (1, 0),
            (4, 2),
            (3, 0),
            (2, 0),
            (1, 0),
            (5, 2),
            //(2, 2),
        ];
        for i in 0..state.len() {
            assert_eq!(gray.current_forward_state(), state[i]);
            gray.next_column();
        }
        println!("current_column: {:?}", gray.current_column());
        println!("current_forward_state: {:?}", gray.current_forward_state());
        println!("forward_index: {:?}", gray.forward_index);
        println!("backward_index: {:?}", gray.backward_index);
        gray.previous_column();
        println!("current_backward_state: {:?}", gray.current_backward_state());
        println!("current_column: {:?}", gray.current_column());
        println!("forward_index: {:?}", gray.forward_index);
        println!("backward_index: {:?}", gray.backward_index);
        gray.next_column();
        println!("current_column: {:?}", gray.current_column());
        println!("current_forward_state: {:?}", gray.current_forward_state());
        println!("forward_index: {:?}", gray.forward_index);
        println!("backward_index: {:?}", gray.backward_index);
        let column = gray.current_column();
        let state = gray.current_forward_state();
        let forward_index = gray.forward_index.clone();
        let backward_index = gray.backward_index.clone();

        for _ in 0..10 {
            gray.previous_column();
            gray.next_column();
            assert_eq!(gray.current_forward_state(), state);
            assert_eq!(gray.forward_index, forward_index);
            assert_eq!(gray.backward_index, backward_index);
            assert_eq!(gray.current_column(), column);
        }
    }

    #[test]
    fn test_gray_code_forward_backward() {
        let w = 5;
        let mut gray = Gray::new(2*w,w);
        let mut state: Vec<(usize, usize)> = Vec::new();
        let mut column: Vec<Vec<u8>> = Vec::new();
        let length = binomial(w)-1;
        for i in 0..length {
            println!("i: {}", i);
            state.push(gray.current_forward_state());
            column.push(gray.current_column());
            println!("current_state: {:?}", gray.current_forward_state());
            println!("current_column: {:?}", gray.current_column());
            gray.next_column();
        }
        for i in 0..length {
            println!("i: {}", length-i-1);
            gray.previous_column();
            //assert_eq!(gray.current_forward_state(), state[length-i-1]);
            println!("current_state: {:?}", gray.current_forward_state());
            println!("current_column: {:?}", gray.current_column());
            assert_eq!(gray.current_column(), column[length-i-1]);
        }
    }

    #[test]
    fn test_gray_code_next_column() {
        let mut gray = Gray::new(6, 3);
        let _expected_states = [
            (3, 3),
            (1, 1),
            (2, 1),
            (4, 3),
            (2, 1),
            (5, 3),
        ];
        for _ in 0..63 {
            let current_state = gray.current_forward_state();
            println!("current_state: {:?}", current_state);
            println!("current_column: {:?}", gray.current_column());
            gray.next_column();
        }
    }
}