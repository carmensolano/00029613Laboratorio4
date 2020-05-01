# Laboratorio 4
Técnicas de Simulación en Computadoras: Cuarta práctica de laboratorio 

## Contenido

<p align="center"> Método de los Elementos Finitos para el problema de dinámica de fluidos en una dimensión 
con funciones de forma lineales y con peso de Galerkin </p>

### Cambios en classes.h

```cpp
//Se actualiza esta enumeracion para tomar en cuenta los nuevos parametros
enum parameters {ELEMENT_LENGTH,ADJECTIVE_VELOCITY,DYNAMIC_VISCOSITY,DENSITY,EXTERNAL_FORCE};

//Esta enumeracion ya no necesita el valor de Neumann
enum sizes {NODES,ELEMENTS,DIRICHLET};
```

```cpp
class item{
    protected:
        int id;
        float x;
        int node1;
        int node2;
        float value;
    public:
    //Se añaden setters para los atributos
        void setId(int identifier) {
            id = identifier;
        }

        void setX(float x_coord) {
            x = x_coord;
        }

        void setNode1(int node_1) {
            node1 = node_1;
        }

        void setNode2(int node_2) {
            node2 = node_2;
        }

        void setValue(float value_to_assign) {
            value = value_to_assign;
        }
```

```cpp
class mesh{
        //Se modifican los tamaños de los arreglos para adecuarse a las nuevas circumstancias
        float parameters[5];
        int sizes[3];
        //Se elimina el arreglo para las condiciones de Neumann
        node *node_list;
        element *element_list;
        condition *dirichlet_list;
    public:
        //Se actualizan los metodos para reflejar la nueva cantidad de parametros
        //y de arreglos
        void setParameters(float l,float u_bar,float nu,float rho,float f){
            parameters[ELEMENT_LENGTH]=l;
            parameters[ADJECTIVE_VELOCITY]=u_bar;
            parameters[DYNAMIC_VISCOSITY]=nu;
            parameters[DENSITY]=rho;
            parameters[EXTERNAL_FORCE]=f;
        }
   }
```

### Cambios en problem.msh

La primera línea corresponde a respectivamente

```
0.3 1.4 0.5 1000 10
10 9 2 1
Coordinates
```

Adicionalmente, tenemos condiciones para cada variable en cada nodo

```
DirichletU
1	15
10  3
EndDirichletU

DirichletP
10  0
EndDirichletP
```

### Cambios en sel.h

Codificamos nuestras nuevas matrices

```cpp
//Matriz resultante de la aplicacion del MEF al termino convectivo
void createLocalA(Matrix &A,mesh m){
    float u_bar = m.getParameter(ADJECTIVE_VELOCITY);
    A.at(0).at(0) += -u_bar/2;  A.at(0).at(1) += u_bar/2;
    A.at(1).at(0) += -u_bar/2;  A.at(1).at(1) += u_bar/2;
}

//Matriz resultante de la aplicacion del MEF al termino relacionado a la viscosidad
void createLocalB(Matrix &B,mesh m){
    float l = m.getParameter(ELEMENT_LENGTH);
    float nu = m.getParameter(DYNAMIC_VISCOSITY);
    B.at(0).at(0) += nu/l;      B.at(0).at(1) += -nu/l;
    B.at(1).at(0) += -nu/l;     B.at(1).at(1) += nu/l;
}

//Matriz resultante de la aplicacion del MEF al termino de la presion
void createLocalC(Matrix &C,mesh m){
    float rho = m.getParameter(DENSITY);
    C.at(0).at(0) += -1/(2*rho);    C.at(0).at(1) += 1/(2*rho);
    C.at(1).at(0) += -1/(2*rho);    C.at(1).at(1) += 1/(2*rho);
}

//Matriz resultante de la aplicacion del MEF al termino de la divergencia
void createLocalD(Matrix &D,mesh m){
    D.at(0).at(0) += -0.5;  D.at(0).at(1) += 0.5;
    D.at(1).at(0) += -0.5;  D.at(1).at(1) += 0.5;
}
```

Se actualiza la funcion para reflejar la nueva estructura de la matrix local K

```cpp
Matrix createLocalK(int element,mesh &m){
    Matrix K,A,B,C,D;

    zeroes(A,2);
    zeroes(B,2);
    zeroes(C,2);
    zeroes(D,2);
    createLocalA(A,m);
    createLocalB(B,m);
    createLocalC(C,m);
    createLocalD(D,m);

    Vector row1, row2, row3, row4;
}
```

Se construyen las filas de la matriz K de acuerdo a su configuracion interna

![equation](https://latex.codecogs.com/gif.latex\begin{pmatrix}&space;A&plus;B&space;&&space;C\\&space;D&space;&&space;0&space;\end{pmatrix})

```cpp
    row1.push_back(A.at(0).at(0)+B.at(0).at(0)); 
    row1.push_back(A.at(0).at(1)+B.at(0).at(1));
    row1.push_back(C.at(0).at(0));                  
    row1.push_back(C.at(0).at(1));

    row2.push_back(A.at(1).at(0)+B.at(1).at(0)); 
    row2.push_back(A.at(1).at(1)+B.at(1).at(1));
    row2.push_back(C.at(1).at(0)); 
    row2.push_back(C.at(1).at(1));

    row3.push_back(D.at(0).at(0)); 
    row3.push_back(D.at(0).at(1));
    row3.push_back(0); 
    row3.push_back(0);

    row4.push_back(D.at(1).at(0)); 
    row4.push_back(D.at(1).at(1));
    row4.push_back(0); 
    row4.push_back(0);

    K.push_back(row1); 
    K.push_back(row2); 
    K.push_back(row3); 
    K.push_back(row4);
```

Se actualiza la construccion de la b local de acuerdo a su nueva estructura

```cpp
Vector createLocalb(int element,mesh &m){
    Vector b;

    float f = m.getParameter(EXTERNAL_FORCE), l = m.getParameter(ELEMENT_LENGTH);
    
    b.push_back(f*l/2); 
    b.push_back(f*l/2); 
    b.push_back(0); 
    b.push_back(0);

    return b;
}
```
Al momento de ensamblar, para los índices de la presión se toma en cuenta el desplazamiento al estar todas las presiones 
después de todas las velocidades

```cpp
void assemblyK(element e,Matrix localK,Matrix &K,int nnodes){
    int index1 = e.getNode1() - 1;
    int index2 = e.getNode2() - 1;
    int index3 = index1 + nnodes;
    int index4 = index2 + nnodes;
}
```

### Cambios en tools.h
Se preparan localmente arreglos para guardar por separado las condiciones de Dirichlet para la velocidad y las de la presion.

```cpp
condition *dirichlet_u_list;
condition *dirichlet_p_list;
```
Se actualizan los datos de las condiciones de Dirichlet de la presión
Se unen los arreglos locales de las condiciones de Dirichlet
Se corrijen los valores de las condiciones

```cpp
updateConditionNodes(ndirich_p,dirichlet_p_list,nnodes);
joinConditions(m.getDirichlet(),ndirich_u,ndirich_p,dirichlet_u_list,dirichlet_p_list);
correctConditions(ndirich_u+ndirich_p,m.getDirichlet());
```

Se agrega una funcion que actualice el indice de todas las condiciones de Dirichlet
para la presion, tomando en cuenta que si una condicion esta en el nodo p, en la matriz
global correspondera a la posicion p+n, donde n es el total de nodos.

```cpp
void updateConditionNodes(int n,condition * list,int delta){
    for(int i=0;i<n;i++)
        list[i].setNode1(list[i].getNode1()+delta);
}
```

Se agrega una función que una las condiciones de Dirichlet para la velocidad, y las condiciones de Dirichlet para la presión, de forma que sean un solo conjunto de condiciones.

```cpp
void joinConditions(condition *list,int n1,int n2,condition *list1,condition *list2){
    int i;
    for(i=0;i<n1;i++)
        list[i] = list1[i];
    for(int j=0;j<n2;j++){
        list[i] = list2[j];
        i++;
    }
}
```
Se agrega una función que modifique los nodos de las condiciones de Dirichlet, de forma que cada
condición tome en cuenta la aplicación de la condición anterior.

```cpp
void correctConditions(int n,condition *list){
    for(int i=0;i<n-1;i++){
        int pivot = list[i].getNode1();
        for(int j=i;j<n;j++)
            /*Si la condición actual corresponde a un nodo posterior al nodo eliminado por
            aplicar la condición anterior, se debe actualizar su posición.*/
            if(list[j].getNode1()>pivot)
                list[j].setNode1(list[j].getNode1()-1);
    }
}
```

<hr>
<p align="center">Para servirles, <strong>Equipo de Instructores</strong> </p>
