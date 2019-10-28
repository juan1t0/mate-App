#include "dyrray.h"

#define MAX(x, y) (((x) > (y)) ? (x) : (y))

int* new_array(int* ptr_last){
    int * arr = malloc(1*sizeof(int));
    ptr_last = arr;
    return arr;
}

void realloc_array(int * arr, unsigned int newSize, int * ptr_last){
    arr = realloc(arr,newSize * sizeof(int));
    if((ptr_last - arr) > newSize)
        ptr_last = arr + (newSize-1);
}

bool array_is_empty(int *arr, int* ptr_last){
    return (ptr_last < arr);
}

void push_element(int* arr, int element, int *ptr_last){
    unsigned int size = ptr_last - arr;
    realloc_array(arr,size+1,ptr_last);    
    *ptr_last = element;
}

int get_element(int* arr, int index, int *ptr_last){
    unsigned int size = ptr_last - arr;
    if(index >= size) return -1;
    return *(arr + index);
}

void set_element(int* arr, int index, int element, int *ptr_last){
    unsigned int size = ptr_last - arr;
    if(index >= size) return -1;
    *(arr + index) = element;
}

void print_array(int* arr, int *ptr_last){
    unsigned int size = ptr_last - arr;
    printf("[");
    if(!array_is_empty(arr,ptr_last)){
        printf("%f", *(arr));
        for(int i=1; i<size; ++i) {
            printf(", %f",*(arr + i));
        }
    }
    printf("]\n");
}

void delete_element(int* arr, int index, int *ptr_last){
    unsigned int size = ptr_last - arr;
    if (index >= size) return;
    //int *aux = arr + index;
    for(int i = index;i<size-1;++i){
        *(arr + i) = *(arr + i + 1);
    }
    ptr_last--;
}

void Intersec(int *targe,int * last_targe,int * A,int *last_A,int *B,int *last_B){
    int *temp = new_array(last_targe);
    unsigned int size_A = last_A - A;
    unsigned int size_B = last_B - B;
    for(int i=0;i<size_A;++i){
        for(int j=0;j<size_B;++j){
            if(*(A+i) == *(B+j))
                push_element(temp,*(A+i),last_targe);
        }
    }
    targe = temp;
}