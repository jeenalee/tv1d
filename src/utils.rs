/// TODO
pub fn sync_values(anchor_value: usize, values: &mut [&mut usize]) -> () {
    for value in values {
        **value = anchor_value;
    }
}

use core::iter::FromIterator;
use core::{iter,mem,slice};

///Fixed size circular/cyclic/ring buffer
///
///A FIFO (first in, first out) queue.
///It cannot represent an empty buffer.
///
///When constructed, the internal `list` must not be empty, and cannot contain invalid (e.g. uninitialized) elements.
#[derive(Clone,Eq,PartialEq,Hash,Debug)]
pub struct CircularBuffer<T>{
    list: Box<[T]>,
    first: usize,
}

impl<T> CircularBuffer<T>{
    pub fn first(&self) -> usize {
        self.first
    }

    ///Returns the number of elements (before starting to loop around).
    #[inline]
    pub fn len(&self) -> usize{self.list.len()}

    ///Enqueues (push at beginning) the given element at the beginning of the buffer
    ///Dequeues (pop at end) the last element and returns it
    ///This keeps the the buffer length
    pub fn rev_queue(&mut self,mut elem: T) -> T{
        let len = self.len();
        mem::swap(
            unsafe {
                self.list.get_unchecked_mut(self.first)
            }
            ,&mut elem);
        self.first = (self.first + 1) % len;
        elem
    }

    ///Enqueues (push at beginning) the given element at the beginning of the buffer
    ///Dequeues (pop at end) the last element and returns it
    ///This keeps the the buffer length
    pub fn queue(&mut self,mut elem: T) -> T{
	let len = self.len();
	self.first = (self.first + len - 1) % len;
	mem::swap(unsafe{self.list.get_unchecked_mut(self.first)},&mut elem);
	elem
    }

//     ///Sets the offset for the first element, relative to the currently first element
//     ///When `index` is out of range, it loops around
//     pub fn set_first(&mut self,index: usize){
// 	self.first = (index + self.first) % self.len();
//     }

//     ///Returns a reference to the element at the given index
//     ///When `index` is out of range, it loops around
//     pub fn get(&self,index: usize) -> &T{
// 	let len = self.len();
// 	unsafe{self.list.get_unchecked((index + self.first) % len)}
//     }

//     ///Returns a mutable reference to the element at the given index
//     ///When `index` is out of range, it loops around
//     pub fn get_mut(&mut self,index: usize) -> &mut T{
// 	let len = self.len();
// 	unsafe{self.list.get_unchecked_mut((index + self.first) % len)}
//     }

//     ///Swaps the two elements at the given indices `a` and `b`.
//     ///When `a` or `b` are out of range, they loop around
//     pub fn swap_internal(&mut self,a: usize,b: usize){
// 	let len = self.len();
// 	self.list.swap((a + self.first) % len,(b + self.first) % len);
//     }

//     ///Swaps the element at the given index with the specifiied new one.
//     ///When `a` or `b` are out of range, they loop around
//     pub fn swap(&mut self,index: usize,mut elem: T) -> T{
// 	mem::swap(self.get_mut(index),&mut elem);
// 	elem
//     }

//     ///Returns an iterator over the buffer looping around at the end.
//     ///This creates a never ending iterator
//     pub fn iter_circular<'s>(&'s self) -> IterCircular<'s,T>{
// 	self.list.iter().cycle().skip(self.first)
//     }

//     ///Returns an iterator over the buffer without looping around.
//     pub fn iter<'s>(&'s self) -> Iter<'s,T>{
// 	self.iter_circular().take(self.len())
//     }
}

impl<T> From<Vec<T>> for CircularBuffer<T>{
	#[inline]
	fn from(vec: Vec<T>) -> Self{
	    debug_assert!(vec.len() > 0);
            CircularBuffer{
		list: vec.into_boxed_slice(),
		first: 0,
            }
	}
}

pub type Iter<'t,T> = iter::Take<IterCircular<'t,T>>;
pub type IterCircular<'t,T> = iter::Skip<iter::Cycle<slice::Iter<'t,T>>>;

#[test]
fn test_len(){
	let l = CircularBuffer::from(Box::new(['a','b','c','d']) as Box<[char]>);
	assert_eq!(l.len(),4);

	let l = CircularBuffer::from(Box::new(['a','b']) as Box<[char]>);
	assert_eq!(l.len(),2);

	let l = CircularBuffer::from(Box::new(['a']) as Box<[char]>);
	assert_eq!(l.len(),1);
}

#[test]
#[should_panic]
fn test_len_empty(){
	let _ = CircularBuffer::from(Box::new([]) as Box<[char]>);
}

#[test]
fn test_queue(){
	let mut l = CircularBuffer::from(Box::new(['a','b','c','d']) as Box<[char]>);
	assert_eq!(l.first,0);
	assert_eq!(&*l.list,&['a','b','c','d']);

	l.queue('9');
	assert_eq!(l.first,3);
	assert_eq!(&*l.list,&['a','b','c','9']);

	l.queue('8');
	assert_eq!(l.first,2);
	assert_eq!(&*l.list,&['a','b','8','9']);

	l.queue('7');
	assert_eq!(l.first,1);
	assert_eq!(&*l.list,&['a','7','8','9']);

	l.queue('6');
	assert_eq!(l.first,0);
	assert_eq!(&*l.list,&['6','7','8','9']);

	l.queue('5');
	assert_eq!(l.first,3);
	assert_eq!(&*l.list,&['6','7','8','5']);

	l.queue('4');
	assert_eq!(l.first,2);
	assert_eq!(&*l.list,&['6','7','4','5']);
}

#[test]
fn test_set_first(){
	let mut l = CircularBuffer::from(Box::new(['a','b','c','d']) as Box<[char]>);
	l.set_first(0);
	assert_eq!(l.first,0);

	l.set_first(1);
	assert_eq!(l.first,1);

	l.set_first(1);
	assert_eq!(l.first,2);

	l.set_first(1);
	assert_eq!(l.first,3);

	l.set_first(1);
	assert_eq!(l.first,0);

	l.set_first(2);
	assert_eq!(l.first,2);

	l.set_first(4);
	assert_eq!(l.first,2);
}

#[test]
fn test_get(){
	let mut l = CircularBuffer::from(Box::new(['a','b','c','d']) as Box<[char]>);
	l.set_first(0);
	assert_eq!(l.first,0);
	assert_eq!(*l.get(0),'a');
	assert_eq!(*l.get(1),'b');
	assert_eq!(*l.get(2),'c');
	assert_eq!(*l.get(3),'d');

	let mut l = CircularBuffer::from(Box::new(['a','b','c','d']) as Box<[char]>);
	l.set_first(1);
	assert_eq!(l.first,1);
	assert_eq!(*l.get(0),'b');
	assert_eq!(*l.get(1),'c');
	assert_eq!(*l.get(2),'d');
	assert_eq!(*l.get(3),'a');

	let mut l = CircularBuffer::from(Box::new(['a','b','c','d']) as Box<[char]>);
	l.set_first(2);
	assert_eq!(l.first,2);
	assert_eq!(*l.get(0),'c');
	assert_eq!(*l.get(1),'d');
	assert_eq!(*l.get(2),'a');
	assert_eq!(*l.get(3),'b');

	let mut l = CircularBuffer::from(Box::new(['a','b','c','d']) as Box<[char]>);
	l.set_first(3);
	assert_eq!(l.first,3);
	assert_eq!(*l.get(0),'d');
	assert_eq!(*l.get(1),'a');
	assert_eq!(*l.get(2),'b');
	assert_eq!(*l.get(3),'c');

	let mut l = CircularBuffer::from(Box::new(['a','b','c','d']) as Box<[char]>);
	l.set_first(4);
	assert_eq!(l.first,0);
	assert_eq!(*l.get(0),'a');
	assert_eq!(*l.get(1),'b');
	assert_eq!(*l.get(2),'c');
	assert_eq!(*l.get(3),'d');

	let mut l = CircularBuffer::from(Box::new(['a','b','c','d']) as Box<[char]>);
	l.set_first(5);
	assert_eq!(l.first,1);
	assert_eq!(*l.get(0),'b');
	assert_eq!(*l.get(1),'c');
	assert_eq!(*l.get(2),'d');
	assert_eq!(*l.get(3),'a');
}

#[test]
fn test_get_mut(){
	let mut l = CircularBuffer::from(Box::new(['a','b','c','d']) as Box<[char]>);
	l.set_first(0);
	assert_eq!(l.first,0);
	assert_eq!(*l.get_mut(0),'a');
	assert_eq!(*l.get_mut(1),'b');
	assert_eq!(*l.get_mut(2),'c');
	assert_eq!(*l.get_mut(3),'d');

	let mut l = CircularBuffer::from(Box::new(['a','b','c','d']) as Box<[char]>);
	l.set_first(1);
	assert_eq!(l.first,1);
	assert_eq!(*l.get_mut(0),'b');
	assert_eq!(*l.get_mut(1),'c');
	assert_eq!(*l.get_mut(2),'d');
	assert_eq!(*l.get_mut(3),'a');

	let mut l = CircularBuffer::from(Box::new(['a','b','c','d']) as Box<[char]>);
	l.set_first(2);
	assert_eq!(l.first,2);
	assert_eq!(*l.get_mut(0),'c');
	assert_eq!(*l.get_mut(1),'d');
	assert_eq!(*l.get_mut(2),'a');
	assert_eq!(*l.get_mut(3),'b');

	let mut l = CircularBuffer::from(Box::new(['a','b','c','d']) as Box<[char]>);
	l.set_first(3);
	assert_eq!(l.first,3);
	assert_eq!(*l.get_mut(0),'d');
	assert_eq!(*l.get_mut(1),'a');
	assert_eq!(*l.get_mut(2),'b');
	assert_eq!(*l.get_mut(3),'c');

	let mut l = CircularBuffer::from(Box::new(['a','b','c','d']) as Box<[char]>);
	l.set_first(4);
	assert_eq!(l.first,0);
	assert_eq!(*l.get_mut(0),'a');
	assert_eq!(*l.get_mut(1),'b');
	assert_eq!(*l.get_mut(2),'c');
	assert_eq!(*l.get_mut(3),'d');

	let mut l = CircularBuffer::from(Box::new(['a','b','c','d']) as Box<[char]>);
	l.set_first(5);
	assert_eq!(l.first,1);
	assert_eq!(*l.get_mut(0),'b');
	assert_eq!(*l.get_mut(1),'c');
	assert_eq!(*l.get_mut(2),'d');
	assert_eq!(*l.get_mut(3),'a');
}

#[test]
fn test_swap(){
	let mut l = CircularBuffer::from(Box::new(['a','b','c','d']) as Box<[char]>);
	assert_eq!(&*l.list,&['a','b','c','d']);

	l.swap(0,'0');
	assert_eq!(&*l.list,&['0','b','c','d']);

	l.swap(1,'1');
	assert_eq!(&*l.list,&['0','1','c','d']);

	l.swap(2,'2');
	assert_eq!(&*l.list,&['0','1','2','d']);

	l.swap(3,'3');
	assert_eq!(&*l.list,&['0','1','2','3']);

	l.swap(4,'4');
	assert_eq!(&*l.list,&['4','1','2','3']);

	l.swap(5,'5');
	assert_eq!(&*l.list,&['4','5','2','3']);

	let mut l = CircularBuffer::from(Box::new(['a','b','c','d']) as Box<[char]>);
	l.set_first(1);
	assert_eq!(&*l.list,&['a','b','c','d']);

	l.swap(0,'0');
	assert_eq!(&*l.list,&['a','0','c','d']);

	l.swap(1,'1');
	assert_eq!(&*l.list,&['a','0','1','d']);

	l.swap(2,'2');
	assert_eq!(&*l.list,&['a','0','1','2']);

	l.swap(3,'3');
	assert_eq!(&*l.list,&['3','0','1','2']);

	l.swap(4,'4');
	assert_eq!(&*l.list,&['3','4','1','2']);

	l.swap(5,'5');
	assert_eq!(&*l.list,&['3','4','5','2']);
}

#[test]
fn test_swap_internal(){
	let mut l = CircularBuffer::from(Box::new(['a','b','c','d']) as Box<[char]>);
	assert_eq!(&*l.list,&['a','b','c','d']);

	l.swap_internal(0,3);
	assert_eq!(&*l.list,&['d','b','c','a']);

	l.swap_internal(3,0);
	assert_eq!(&*l.list,&['a','b','c','d']);

	l.swap_internal(1,2);
	assert_eq!(&*l.list,&['a','c','b','d']);

	l.swap_internal(2,1);
	assert_eq!(&*l.list,&['a','b','c','d']);

	l.swap_internal(0,5);
	assert_eq!(&*l.list,&['b','a','c','d']);

	l.swap_internal(5,0);
	assert_eq!(&*l.list,&['a','b','c','d']);
}

#[test]
fn test_iter(){
	let l = unsafe{CircularBuffer::from_raw_parts(Box::new(['a','b','c']) as Box<[char]>,0)};
	let mut i = l.iter();

	assert_eq!(*i.next().unwrap(),'a');
	assert_eq!(*i.next().unwrap(),'b');
	assert_eq!(*i.next().unwrap(),'c');
	assert!(i.next().is_none());

	let l = unsafe{CircularBuffer::from_raw_parts(Box::new(['a','b','c']) as Box<[char]>,1)};
	let mut i = l.iter();

	assert_eq!(*i.next().unwrap(),'b');
	assert_eq!(*i.next().unwrap(),'c');
	assert_eq!(*i.next().unwrap(),'a');
	assert!(i.next().is_none());

	let l = unsafe{CircularBuffer::from_raw_parts(Box::new(['a','b','c']) as Box<[char]>,2)};
	let mut i = l.iter();

	assert_eq!(*i.next().unwrap(),'c');
	assert_eq!(*i.next().unwrap(),'a');
	assert_eq!(*i.next().unwrap(),'b');
	assert!(i.next().is_none());
}

#[test]
fn test_iter_circular(){
	let l = unsafe{CircularBuffer::from_raw_parts(Box::new(['a','b','c']) as Box<[char]>,0)};
	let mut i = l.iter_circular();

	assert_eq!(*i.next().unwrap(),'a');
	assert_eq!(*i.next().unwrap(),'b');
	assert_eq!(*i.next().unwrap(),'c');
	assert_eq!(*i.next().unwrap(),'a');
	assert_eq!(*i.next().unwrap(),'b');
	assert_eq!(*i.next().unwrap(),'c');
	assert_eq!(*i.next().unwrap(),'a');
	assert_eq!(*i.next().unwrap(),'b');
	assert_eq!(*i.next().unwrap(),'c');

	let l = unsafe{CircularBuffer::from_raw_parts(Box::new(['a','b','c']) as Box<[char]>,1)};
	let mut i = l.iter_circular();

	assert_eq!(*i.next().unwrap(),'b');
	assert_eq!(*i.next().unwrap(),'c');
	assert_eq!(*i.next().unwrap(),'a');
	assert_eq!(*i.next().unwrap(),'b');
	assert_eq!(*i.next().unwrap(),'c');
	assert_eq!(*i.next().unwrap(),'a');
	assert_eq!(*i.next().unwrap(),'b');
	assert_eq!(*i.next().unwrap(),'c');
	assert_eq!(*i.next().unwrap(),'a');

	let l = unsafe{CircularBuffer::from_raw_parts(Box::new(['a','b','c']) as Box<[char]>,2)};
	let mut i = l.iter_circular();

	assert_eq!(*i.next().unwrap(),'c');
	assert_eq!(*i.next().unwrap(),'a');
	assert_eq!(*i.next().unwrap(),'b');
	assert_eq!(*i.next().unwrap(),'c');
	assert_eq!(*i.next().unwrap(),'a');
	assert_eq!(*i.next().unwrap(),'b');
	assert_eq!(*i.next().unwrap(),'c');
	assert_eq!(*i.next().unwrap(),'a');
	assert_eq!(*i.next().unwrap(),'b');
}
