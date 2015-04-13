//#pragma once
///*	queue.h
//
//	Header file for queue implementation
//
//	by: Steven Skiena
//*/
//
///*
//Copyright 2003 by Steven S. Skiena; all rights reserved. 
//
//Permission is granted for use in non-commerical applications
//provided this copyright notice remains intact and unchanged.
//
//This program appears in my book:
//
//"Programming Challenges: The Programming Contest Training Manual"
//by Steven Skiena and Miguel Revilla, Springer-Verlag, New York 2003.
//
//See our website www.programming-challenges.com for additional information.
//
//This book can be ordered from Amazon.com at
//
//http://www.amazon.com/exec/obidos/ASIN/0387001638/thealgorithmrepo/
//
//*/
//

#define TYPE int

class MyQueue {
	TYPE* q;		/* body of queue */
	int first;                      /* position of first element */
	int last;                       /* position of last element */
	int count;                      /* number of queue elements */
	int size_;

public:
	void init(int size)
	{
		this->size_ = size;
		this->q = new TYPE [this->size_+1];
		this->first = 0;
		this->last = this->size_-1;
		this->count = 0;
	}

	void push(TYPE x)
	{
		if (this->count >= this->size_)
		{
			//printf("Warning: queue overflow enqueue x=%d\n",x);
		}
		else {
			this->last = (this->last+1) % this->size_;
			this->q[ this->last ] = x;    
			this->count = this->count + 1;
		}
	}

	int pop()
	{
		TYPE x;

		if (this->count <= 0) 
		{//printf("Warning: empty queue dequeue.\n");
		}
		else {
			x = this->q[ this->first ];
			this->first = (this->first+1) % this->size_;
			this->count = this->count - 1;
		}

		return(x);
	}
	int front()
	{
		TYPE x;

		if (this->count <= 0) 
		{//printf("Warning: empty queue dequeue.\n");
		}
		else {
			return this->q[ this->first ];
		}

	}

	int empty()
	{
		if (this->count <= 0) 
			return true;
		else 
			return false;
	}
	int size()
	{
		return this->size_;
	}
	~MyQueue()
	{
		delete[] q;
	}
} ;
