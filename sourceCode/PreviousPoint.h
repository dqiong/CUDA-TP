#ifndef __PREVIOUSPOINT_H__
#define __PREVIOUSPOINT_H__

#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <device_functions.h>
#include <thrust/host_vector.h>
#include <cuda.h>

#include <iostream>
#include <math.h>
using namespace std;

//record the previous point if they are on the same diagonal.
struct Point{
	int x,y;
};

struct Node{
	Node* ch[2];
	int r;
	int xx,yy;
	double v;
	double eps = 0.05;
   __host__ __device__ bool operator<(Node& rhs)const{
		return r < rhs.r;
	}
	__host__ __device__ int cmp(double x)const{
		if(fabs(x-v)<eps)return -1;
		return x<v ?0:1;
	}
	__host__ __device__ void Rotate(Node* &o,int d)
	{
		Node* k = o ->ch[d^1];
		o->ch[d^1] = k->ch[d];
		k->ch[d] = o;
		o = k;
	}
	__host__ __device__ Node* Find(Node* o,double x)
	{
		while(o != NULL){
			int d = o->cmp(x);
			if(d == -1)return o;
			else o = o->ch[d];
		}
		return NULL;
	}
	__host__ __device__  void Insert(Node* &o,double x,int xx,int yy)
	{
			if(o == NULL){
				o = new Node();
				o ->ch[0] = o->ch[1] = NULL;
				o->xx = xx,o->yy = yy;
				o->v = x;
				o->r = rand();
			}else{
				int d = o->cmp(x);
				Insert(o->ch[d],x,xx,yy);
				if(o->ch[d] > o) Rotate(o,d^1);
			}
	}
	__host__ __device__ void gao(Node* o)
	{
		if (o->ch[0] != NULL) {
			gao(o->ch[0]);
			delete o->ch[0];
			o->ch[0] = NULL;
		}
		if (o->ch[1] != NULL){
			gao(o->ch[1]);
			delete o->ch[1];
			o->ch[1] = NULL;
		}
	}
	__host__ __device__ void Delete(Node* root)
	{
		if (root == NULL)return;
		gao(root);
		//delete root;
		//root = NULL;
	}

};
#endif