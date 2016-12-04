#include "SimpleObject.h"
#include <vector>
#include <cmath>
#include <iostream>
using namespace std;
#define n1 500
#define n2 1000
#define n3 50000000
namespace SimpleOBJ
{
    
    int v_p_adj[150000][n1];
    int v_e_adj[150000][n2];
    int edge_id[n3][3];
    CSimpleObject::CSimpleObject(void)
    {
        m_nVertices = -1;
        m_nTriangles = -1;
        m_pTriangleList = NULL;
        m_pVertexList = NULL;
    }
    CSimpleObject::~CSimpleObject(void)
    {
        Destroy();
    }
    
    void CSimpleObject::Destroy()
    {
        if(m_pTriangleList)
            delete []m_pTriangleList;
        if(m_pVertexList)
            delete []m_pVertexList;
        
        m_nVertices = -1;
        m_nTriangles = -1;
        m_pVertexList = NULL;
        m_pTriangleList = NULL;
    }
    
    bool CSimpleObject::LoadFromObj(const char* fn)
    {
        FILE* fp = fopen(fn,"r");
        if(fp==NULL)
        {
            printf("Error: Loading %s failed.\n",fn);
            return false;
        }
        else
        {
            if(Parse(fp))
            {
                printf("Loading from %s successfully.\n",fn);
                printf("Vertex Number = %d\n",m_nVertices);
                printf("Triangle Number = %d\n",m_nTriangles);
                fclose(fp);
                return true;
            }
            else
            {
                fclose(fp);
                return false;
            }
        }
    }
    
    int CSimpleObject::matrixInversion(double a[4][4], int n)
    {
        int *is = new int[n];
        int *js = new int[n];
        int i,j,k;
        double d,p;
        for ( k = 0; k < n; k++)
        {
            d = 0.0;
            for (i=k; i<=n-1; i++)
                for (j=k; j<=n-1; j++)
                {
                    p=fabs(a[i][j]);
                    if (p>d) { d=p; is[k]=i; js[k]=j;}
                }
            if ( 0.0 == d )
            {
                free(is); free(js); //printf("err**not inv\n");
                return(0);
            }
            if (is[k]!=k)
                for (j=0; j<=n-1; j++)
                {
                    p=a[k][j];
                    a[k][j]=a[is[k]][j];
                    a[is[k]][j]=p;
                }
            if (js[k]!=k)
                for (i=0; i<=n-1; i++)
                {
                    p=a[i][k];
                    a[i][k]=a[i][js[k]];
                    a[i][js[k]]=p;
                }
            a[k][k] = 1.0/a[k][k];
            for (j=0; j<=n-1; j++)
                if (j!=k)
                {
                    a[k][j] *= a[k][k];
                }
            for (i=0; i<=n-1; i++)
                if (i!=k)
                    for (j=0; j<=n-1; j++)
                        if (j!=k)
                        {
                            a[i][j] -= a[i][k]*a[k][j];
                        }
            for (i=0; i<=n-1; i++)
                if (i!=k)
                {
                    a[i][k] = -a[i][k]*a[k][k];
                }
        }
        for ( k = n-1; k >= 0; k--)
        {
            if (js[k]!=k)
                for (j=0; j<=n-1; j++)
                {
                    p = a[k][j];
                    a[k][j] = a[js[k]][j];
                    a[js[k]][j]=p;
                }
            if (is[k]!=k)
                for (i=0; i<=n-1; i++)
                { 
                    p = a[i][k];
                    a[i][k]=a[i][is[k]];
                    a[i][is[k]] = p;
                }
        }
        free(is); free(js);
        return(1);
    }
        
    Qmetric CSimpleObject::compute_a_Q(int v1)
    {
        Qmetric Q;
        Vec3f* vlist =  m_pVertexList;
        Array<int,3>* tlist = m_pTriangleList;
        int change1, change2;
        for (int k = 0; k < 4; k++)
            for (int l = 0; l < 4; l++)
                Q.value[k][l] = 0;
        Vec3f v1_value = vlist[v1];
        for (int i = 1; i < v_p_adj[v1][0] + 1; i++)
        {
            if ((tlist[v_p_adj[v1][i]][0] == tlist[v_p_adj[v1][i]][1])
                || (tlist[v_p_adj[v1][i]][0] == tlist[v_p_adj[v1][i]][2])
                || (tlist[v_p_adj[v1][i]][1] == tlist[v_p_adj[v1][i]][2]))
                continue;
            for (int j = 0; j < 3; j++)
            {
                if (v1 == tlist[v_p_adj[v1][i]][j])
                {
                    change1 = tlist[v_p_adj[v1][i]][(j+1)%3];
                    change2 = tlist[v_p_adj[v1][i]][(j+2)%3];
                    double a[4], x1 = v1_value[0], y1 = v1_value[1], z1 = v1_value[2],
                    x2 = vlist[change1][0], y2 = vlist[change1][1], z2 = vlist[change1][2],
                    x3 = vlist[change2][0], y3 = vlist[change2][1], z3 = vlist[change2][2];
                    
                    a[0] = (y2-y1)*(z3-z1) - (z2-z1)*(y3-y1);
                    a[1] = (x3-x1)*(z2-z1) - (x2-x1)*(z3-z1);
                    a[2] = (x2-x1)*(y3-y1) - (x3-x1)*(y2-y1);
                    a[3] = x1*((y3-y1)*(z2-z1)-(y2-y1)*(z3-z1))
                    + y1*((x2-x1)*(z3-z1)-(x3-x1)*(z2-z1))
                    + z1*((x3-x1)*(y2-y1)-(y3-y1)*(x2-x1));
                    double temp = sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
                    if (temp == 0)
                    {
                        cout << "Error! temp == 0" << endl;
                        continue;
                    }
                    a[0] /= temp;
                    a[1] /= temp;
                    a[2] /= temp;
                    a[3] /= temp;
                    for (int k = 0; k < 4; k++)
                        for (int l = 0; l < 4; l++)
                            Q.value[k][l] += a[k] * a[l];
                    break;
                }
                
            }
        }
        return Q;

    }
    
    
    void CSimpleObject::contract()
    {
        /*
         Each element of tlist is a triangle. Each element of vlist is the position of a vertex.
         The three vertices of the i-th phase in tlist is vlist[tlist[i][0]], vlist[tlist[i][1]], vlist[tlist[i][2]].
         */
        Vec3f* vlist =  m_pVertexList;
        Array<int,3>* tlist = m_pTriangleList;
        //Only choose those existing edges
        set<Edge, EdgeSortCriterion>::iterator it;
        it = all_edges.begin();
        while (edge_id[(*it).ID][0] != 1)
        {
            all_edges.erase(it);
            if (all_edges.size() < 1)
            {
                cout << "Heap is empty" << endl;
                exit(0);
            }
            it = all_edges.begin();
        }
        int erase_id = (*it).ID;
        double erase_cost = (*it).cost;
        all_edges.erase(it);
        if (all_edges.size() < 1)
        {
            cout << "Heap is empty" << endl;
            exit(0);
        }
        m_nowTriangles -= 2;
        int erase_v1 = edge_id[erase_id][1], erase_v2 = edge_id[erase_id][2];
        Vec3f target_v = opt_target(erase_v1, erase_v2);
        Vec3f v1_value = vlist[erase_v1], v2_value = vlist[erase_v2];
        m_pVertexList[erase_v1] = target_v;
        m_pVertexList[erase_v2] = target_v;
        v_e_adj[erase_v2][1] = 0;
        for (int i = 1; i < v_p_adj[erase_v2][0] + 1; i++)//是否要做防重复判断？需要
        {
            int flag = 0;
            for (int j = 1; j < v_p_adj[erase_v1][0] + 1; j++)
            {
                if (v_p_adj[erase_v2][i] == v_p_adj[erase_v1][j])
                {
                    flag = 1;
                    break;
                }
            }
            if (flag == 0)
            {
                for (int j = 1; j < v_p_adj[erase_v1][0] + 1; j++)
                {
                    
                    if (((tlist[v_p_adj[erase_v1][j]][0] == tlist[v_p_adj[erase_v1][j]][1])
                         || (tlist[v_p_adj[erase_v1][j]][0] == tlist[v_p_adj[erase_v1][j]][2])
                         || (tlist[v_p_adj[erase_v1][j]][1] == tlist[v_p_adj[erase_v1][j]][2])))
                    {
                        if (~((tlist[v_p_adj[erase_v2][i]][0] == tlist[v_p_adj[erase_v2][i]][1])
                              || (tlist[v_p_adj[erase_v2][i]][0] == tlist[v_p_adj[erase_v2][i]][2])
                              || (tlist[v_p_adj[erase_v2][i]][1] == tlist[v_p_adj[erase_v2][i]][2])))
                        {
                            v_p_adj[erase_v1][j] = v_p_adj[erase_v2][i];
                            flag = 1;
                            break;
                        }
                    }
                }
            }
            if (flag == 0)
            {
                if (~((tlist[v_p_adj[erase_v2][i]][0] == tlist[v_p_adj[erase_v2][i]][1])
                      || (tlist[v_p_adj[erase_v2][i]][0] == tlist[v_p_adj[erase_v2][i]][2])
                      || (tlist[v_p_adj[erase_v2][i]][1] == tlist[v_p_adj[erase_v2][i]][2])))
                {
                    v_p_adj[erase_v1][v_p_adj[erase_v1][0] + 1] = v_p_adj[erase_v2][i];
                    v_p_adj[erase_v1][0] ++;
                    if (v_p_adj[erase_v1][0] > max_face_num)
                        max_face_num = v_p_adj[erase_v1][0];
                }
                if (v_p_adj[erase_v1][0] > n1 - 2)
                {
                    cout << "Too many phases!" << endl;
                    exit(0);
                }
            }
            for (int j = 0; j < 3; j++)
                if (m_pTriangleList[v_p_adj[erase_v2][i]][j] == erase_v2)
                    m_pTriangleList[v_p_adj[erase_v2][i]][j] = erase_v1;
        }
        
        
        for (int i = 2; i < v_e_adj[erase_v1][0] + 2; i++)
        {
            if (edge_id[v_e_adj[erase_v1][i]][0] == 1)//当该边有效
            {
                edge_id[v_e_adj[erase_v1][i]][0] = 0;//将该边设为无效
                int another_v = (edge_id[v_e_adj[erase_v1][i]][1] != erase_v1) ? edge_id[v_e_adj[erase_v1][i]][1] : edge_id[v_e_adj[erase_v1][i]][2];//该边的另一端点
                if (edge_id[v_e_adj[erase_v1][i]][1] != erase_v2 && edge_id[v_e_adj[erase_v1][i]][2] != erase_v2
                    && v_e_adj[another_v][1] == 1)//当另一端点不是v2，另一端点有效
                {
                    edge_id[id_num][0] = 1;
                    edge_id[id_num][1] = edge_id[v_e_adj[erase_v1][i]][1];
                    edge_id[id_num][2] = edge_id[v_e_adj[erase_v1][i]][2];
                    int temp = v_e_adj[erase_v1][i];
                    v_e_adj[erase_v1][i] = id_num;
                    for (int j = 2; j < v_e_adj[another_v][0] + 2; j++)//修改另一端点的该边id
                    {
                        if (v_e_adj[another_v][j] == temp || edge_id[v_e_adj[another_v][j]][0] == 0)
                        {
                            v_e_adj[another_v][j] = id_num;
                            break;
                        }
                    }
                    Edge add_edge;
                    add_edge.ID = id_num;
                    add_edge.cost = comp_err(erase_v1, another_v);
                    all_edges.insert(add_edge);
                    id_num ++;
                    
                    if (id_num >= n3)
                    {
                        cout << "Too many ids!" << endl;
                        exit(0);
                    }
                }
            }
        }
        
        for (int i = 2; i < v_e_adj[erase_v2][0] + 2; i++)
        {
            if (edge_id[v_e_adj[erase_v2][i]][0] == 1)//当该边有效
            {
                edge_id[v_e_adj[erase_v2][i]][0] = 0;//将该边设为无效
                int another_v = (edge_id[v_e_adj[erase_v2][i]][1] != erase_v2) ? edge_id[v_e_adj[erase_v2][i]][1] : edge_id[v_e_adj[erase_v2][i]][2];
                if (edge_id[v_e_adj[erase_v2][i]][1] != erase_v1 && edge_id[v_e_adj[erase_v2][i]][2] != erase_v1
                    && v_e_adj[another_v][1] == 1)//当另一端点不是v1，另一端点有效
                {
                    int temp = v_e_adj[erase_v2][i];
                    if (edge_id[v_e_adj[erase_v2][i]][1] == erase_v2)
                        edge_id[v_e_adj[erase_v2][i]][1] = erase_v1;
                    if (edge_id[v_e_adj[erase_v2][i]][2] == erase_v2)
                        edge_id[v_e_adj[erase_v2][i]][2] = erase_v1;
                    int flag = 0;
                    
                    for (int j = 2; j < v_e_adj[erase_v1][0] + 2; j++)//在v1的现有边中找是否已经与another_v相连
                    {
                        if ((edge_id[v_e_adj[erase_v1][j]][1] == another_v ||
                             edge_id[v_e_adj[erase_v1][j]][2] == another_v) && edge_id[v_e_adj[erase_v1][j]][0] == 1)
                        {
                            flag = 2;
                            break;
                        }
                    }
                    if (flag == 2)//v1若已经和another_v连了有效边
                    {
                        continue;
                    }
                    edge_id[id_num][0] = 1;
                    edge_id[id_num][1] = edge_id[v_e_adj[erase_v2][i]][1];
                    edge_id[id_num][2] = edge_id[v_e_adj[erase_v2][i]][2];
                    for (int j = 2; j < v_e_adj[erase_v1][v_e_adj[erase_v1][0] + 2]; j++)//在v1的现有边中找能被替代的
                    {
                        if (edge_id[v_e_adj[erase_v1][j]][0] == 0)
                        {
                            v_e_adj[erase_v1][j] = id_num;
                            flag = 1;
                            break;
                        }
                    }
                    if (flag == 0) //若v1现有边全部有效
                    {
                        v_e_adj[erase_v1][v_e_adj[erase_v1][0] + 2] = id_num;
                        v_e_adj[erase_v1][0] ++;
                        if (v_e_adj[erase_v1][0] > max_edge_num)
                            max_edge_num = v_e_adj[erase_v1][0];
                    }
                    if (v_e_adj[erase_v1][0] > n2 - 2)
                    {
                        cout << "Too many edges!" << endl;
                        exit(0);
                    }
                    for (int j = 2; j < v_e_adj[another_v][0] + 2; j++)
                    {
                        if (v_e_adj[another_v][j] == temp || edge_id[v_e_adj[another_v][j]][0] == 0)
                        {
                            v_e_adj[another_v][j] = id_num;
                            break;
                        }
                    }
                    
                    Edge add_edge;
                    add_edge.ID = id_num;
                    add_edge.cost = comp_err(erase_v1, another_v);
                    all_edges.insert(add_edge);
                    id_num ++;
                    
                    if (id_num >= n3)
                    {
                        cout << "Too many ids!" << endl;
                        exit(0);
                    }
                }
            }
        }
        for (int i = 2; i < v_e_adj[erase_v1][0] + 2; i++)
        {
            if (edge_id[v_e_adj[erase_v1][i]][0] == 1)//当该边有效
            {
                int another_v = (edge_id[v_e_adj[erase_v1][i]][1] != erase_v1) ? edge_id[v_e_adj[erase_v1][i]][1] : edge_id[v_e_adj[erase_v1][i]][2];
                if (edge_id[v_e_adj[erase_v1][i]][1] != erase_v2 && edge_id[v_e_adj[erase_v1][i]][2] != erase_v2
                    && v_e_adj[another_v][1] == 1)//当另一端点不是v2，另一端点有效
                {
                    for (int k = 2; k < v_e_adj[another_v][0] + 2; k++)
                    {
                        if (edge_id[v_e_adj[another_v][k]][0] == 1)//当该边有效
                        {
                            
                            int another_v2 = (edge_id[v_e_adj[another_v][k]][1] != another_v) ? edge_id[v_e_adj[another_v][k]][1] : edge_id[v_e_adj[another_v][k]][2];
                            if (another_v2 != erase_v1 && another_v2 != erase_v2 &&
                                v_e_adj[another_v2][1] == 1)//当另一端点不是v1\v2，另一端点有效
                            {
                                edge_id[v_e_adj[another_v][k]][0] = 0;
                                edge_id[id_num][0] = 1;
                                edge_id[id_num][1] = edge_id[v_e_adj[another_v][k]][1];
                                edge_id[id_num][2] = edge_id[v_e_adj[another_v][k]][2];
                                int temp = v_e_adj[another_v][k];
                                v_e_adj[another_v][k] = id_num;
                                for (int j = 2; j < v_e_adj[another_v2][0] + 2; j++)//修改另一端点的该边id
                                {
                                    if (v_e_adj[another_v2][j] == temp || edge_id[v_e_adj[another_v2][j]][0] == 0)
                                    {
                                        v_e_adj[another_v2][j] = id_num;
                                        break;
                                    }
                                }
                                Edge add_edge;
                                add_edge.ID = id_num;
                                add_edge.cost = comp_err(another_v2, another_v);
                                all_edges.insert(add_edge);
                                id_num ++;
                                
                                if (id_num >= n3)
                                {
                                    cout << "Too many ids!" << endl;
                                    exit(0);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    
    double CSimpleObject::comp_err(int i, int j)
    {
        Vec3f v1 = m_pVertexList[i];
        Vec3f v2 = m_pVertexList[j];
        double Q[4][4];
        double error = 0;
        Vec3f v = opt_target(i, j);
        Qmetric Q1, Q2;
        Q1 = compute_a_Q(i);
        Q2 = compute_a_Q(j);
        for (int m = 0; m < 4; m++)
            for (int n = 0; n < 4; n++)
                Q[m][n] = Q1.value[m][n] + Q2.value[m][n];
        error = Q[0][0]*v.x*v.x + 2*Q[0][1]*v.x*v.y + 2*Q[0][2]*v.x*v.z + 2*Q[0][3]*v.x
        + Q[1][1]*v.y*v.y + 2*Q[1][2]*v.y*v.z + 2*Q[1][3]*v.y
        + Q[2][2]*v.z*v.z + 2*Q[2][3]*v.z
        + Q[3][3];
        return error;
    }
    
    Vec3f CSimpleObject::opt_target(int i, int j)
    {
        Vec3f v;
        Vec3f v1 = m_pVertexList[i];
        Vec3f v2 = m_pVertexList[j];
        double Q[4][4];
        double error = 0;
        Qmetric Q1, Q2;
        Q1 = compute_a_Q(i);
        Q2 = compute_a_Q(j);
        for (int m = 0; m < 3; m++)
            for (int n = 0; n < 4; n++)
                Q[m][n] = Q1.value[m][n] + Q2.value[m][n];
        
        for(int j = 0; j < 3; j++){
            Q[3][j] = 0;
        }
        Q[3][3] = 1;
        int flag = matrixInversion(Q, 4);
        if(flag==1){
            Vec3f temp(Q[0][3], Q[1][3], Q[2][3]);
            v = temp;
        }
        else{
            //TODO shoose optimal on the segment
            
            v = (v1 + v2) * 0.5;
        }

        return v;
    }
    
    
    bool CSimpleObject::Parse(FILE* fp)
    {
        
        char buf[256];
        int nVertices,nTriangles;
        std::vector<Vec3f>          vecVertices;
        std::vector<Array<int,3> >  vecTriangles;
        
        nVertices = 0;
        nTriangles = 0;
        vecVertices.clear();
        vecTriangles.clear();
        int lineNumber = 0;
        
        while(fscanf(fp, "%s", buf) != EOF)
        {
            lineNumber ++;
            switch(buf[0])
            {
                case '#':				/* comment */
                    /* eat up rest of line */
                    fgets(buf, sizeof(buf), fp);
                    break;
                case 'v':				/* v, vn, vt */
                    switch(buf[1])
                {
                    case '\0':			    /* vertex */
                    {
                        Vec3f vP;
                        if(fscanf(fp, "%f %f %f",
                                  &vP.x,
                                  &vP.y,
                                  &vP.z)==3)
                        {
                            nVertices++;
                            vecVertices.push_back(vP);
                        }
                        else
                        {
                            printf("Error: Wrong Number of Values(Should be 3). at Line %d\n",lineNumber);
                            return false;
                        }
                    }
                        break;
                    default:
                        /* eat up rest of line */
                        fgets(buf, sizeof(buf), fp);
                        break;
                }
                    break;
                    
                case 'f':				/* face */
                {
                    int v,n,t;
                    Array<int,3> vIndices;
                    if(fscanf(fp, "%s", buf)!=1)
                    {
                        printf("Error: Wrong Face at Line %d\n",lineNumber);
                        return false;
                    }
                    
                    /* can be one of %d, %d//%d, %d/%d, %d/%d/%d %d//%d */
                    if (strstr(buf, "//"))
                    {
                        /* v//n */
                        if( sscanf(buf, "%d//%d", &vIndices[0],&n) ==2  &&
                           fscanf(fp, "%d//%d", &vIndices[1], &n) ==2  &&
                           fscanf(fp, "%d//%d", &vIndices[2], &n) ==2)
                        {
                            nTriangles++;
                            vecTriangles.push_back(vIndices);
                        }
                        else
                        {
                            printf("Error: Wrong Face at Line %d\n",lineNumber);
                            return false;
                        }
                        
                    }
                    else if (sscanf(buf, "%d/%d/%d", &v, &t, &n) == 3)
                    {
                        /* v/t/n */
                        vIndices[0] = v;
                        if( fscanf(fp, "%d/%d/%d", &vIndices[1], &t, &n) ==3 &&
                           fscanf(fp, "%d/%d/%d", &vIndices[2], &t, &n) ==3 )
                        {
                            nTriangles++;
                            vecTriangles.push_back(vIndices);
                        }
                        else
                        {
                            printf("Error: Wrong Face at Line %d\n",lineNumber);
                            return false;
                        }
                    }
                    else if (sscanf(buf, "%d/%d", &v, &t) == 2)
                    {
                        /* v/t */
                        vIndices[0] = v;
                        if( fscanf(fp, "%d/%d", &vIndices[1], &t) ==2 &&
                           fscanf(fp, "%d/%d", &vIndices[2], &t) ==2 )
                        {
                            nTriangles++;
                            vecTriangles.push_back(vIndices);
                        }
                        else
                        {
                            printf("Error: Wrong Face at Line %d\n",lineNumber);
                            return false;
                        }
                    }
                    else
                    {
                        /* v */
                        if( sscanf(buf, "%d", &vIndices[0]) ==1 &&
                           fscanf(fp, "%d", &vIndices[1])  ==1 &&
                           fscanf(fp, "%d", &vIndices[2])  ==1 )
                        {
                            nTriangles++;
                            vecTriangles.push_back(vIndices);
                        }
                        else
                        {
                            printf("Error: Wrong Face at Line %d\n",lineNumber);
                            return false;
                        }
                    }
                    
                }
                    
                    break;
                    
                default:
                    /* eat up rest of line */
                    fgets(buf, sizeof(buf), fp);
                    break;
            }
        }
        
        if(CheckParse(vecVertices.size(),vecTriangles))
        {
            Destroy();
            
            m_nVertices = vecVertices.size();
            m_nTriangles = vecTriangles.size();
            m_pVertexList = new Vec3f[m_nVertices];
            m_pTriangleList = new Array<int,3> [m_nTriangles];
            m_nowTriangles = m_nTriangles;
            for(int i=0;i<m_nVertices;i++)
                m_pVertexList[i] = vecVertices[i];
            for(int i=0;i<m_nTriangles;i++)
            {
                m_pTriangleList[i][0] = vecTriangles[i][0] - 1;
                m_pTriangleList[i][1] = vecTriangles[i][1] - 1;
                m_pTriangleList[i][2] = vecTriangles[i][2] - 1;
            }
            init_adjlist();
            int t;
            t = 1;
            return true;
        }
        else
            return false;
    }
    
    void CSimpleObject::init_adjlist()
    {
        Vec3f* vlist =  m_pVertexList;
        Array<int,3>* tlist = m_pTriangleList;
        for (int i = 0; i < m_nTriangles; i++)
        {
            v_p_adj[tlist[i][0]][v_p_adj[tlist[i][0]][0] + 1] = i;
            v_p_adj[tlist[i][0]][0] ++;
            v_p_adj[tlist[i][1]][v_p_adj[tlist[i][1]][0] + 1] = i;
            v_p_adj[tlist[i][1]][0] ++;
            v_p_adj[tlist[i][2]][v_p_adj[tlist[i][2]][0] + 1] = i;
            v_p_adj[tlist[i][2]][0] ++;
            for (int j = 0; j < 3; j++)
            {
                int flag = 0;
                for (int k = 2; k < v_e_adj[tlist[i][j]][0] + 2; k++)
                {
                    if (edge_id[v_e_adj[tlist[i][j]][k]][1] == tlist[i][(j+1)%3] || edge_id[v_e_adj[tlist[i][j]][k]][2] == tlist[i][(j+1)%3])
                    {
                        flag = 1;
                        break;
                    }
                }
                if (flag == 0)
                {
                    v_e_adj[tlist[i][j]][v_e_adj[tlist[i][j]][0] + 2] = id_num;
                    v_e_adj[tlist[i][(j+1)%3]][v_e_adj[tlist[i][(j+1)%3]][0] + 2] = id_num;
                    v_e_adj[tlist[i][j]][0] ++;
                    v_e_adj[tlist[i][(j+1)%3]][0] ++;
                    v_e_adj[tlist[i][(j+1)%3]][1]  = 1;
                    v_e_adj[tlist[i][j]][1] = 1;
                    edge_id[id_num][0] = 1;
                    edge_id[id_num][1] = tlist[i][j];
                    edge_id[id_num][2] = tlist[i][(j+1)%3];
                    id_num ++;
                }
            }
        }
        for (int i = 0; i < id_num; i++)
        {
            Edge e;
            e.ID = i;
            e.cost = comp_err(edge_id[i][1], edge_id[i][2]);
            if (e.cost < 0)
            {
                cout << "error! cost less than 0" << endl;
            }
            double test;
            test = comp_err(edge_id[i][1], edge_id[i][2]);
            all_edges.insert(e);
        }
    }
    
    bool CSimpleObject::CheckParse(int nVertices,std::vector<Array<int,3> > & vecTriangles)
    {
        for(int i=0;i<vecTriangles.size();i++)
        {
            Array<int,3> & vIndices = vecTriangles[i];
            for(int j=0;j<vIndices._len;j++)
                if(vIndices[j]<=0||vIndices[j]>nVertices)
                {
                    printf("Error: The vertex index of Face %d has exceeded vertex number %d\n",i,nVertices);
                    return false;
                }
        }
        
        return true;
    }
    
    void CSimpleObject::edit()
    {
        for (int i = 0; i < m_nVertices; i++)//对每一个有效点
        {
            if (v_e_adj[i][1] != 0)
            {
                for (int j = 0; j < v_p_adj[i][0] + 1; j++)
                {
                    int face = v_p_adj[i][j];
                    for (int k = 0; k < 3; k++)
                    {
                        if (m_pTriangleList[face][k] == i)
                        {
                            m_pTriangleList[face][k] = m_newVertex;
                        }
                    }
                }
                m_pVertexList[m_newVertex] = m_pVertexList[i];
                m_newVertex++;
            }
        }
        for (int i = 0; i < m_nTriangles; i++)
        {
            if (~(m_pTriangleList[i][0] == m_pTriangleList[i][1])
                || (m_pTriangleList[i][0] == m_pTriangleList[i][2])
                || (m_pTriangleList[i][1] == m_pTriangleList[i][2]))
                m_newTriangle++;
        }
    }
    
    
    bool CSimpleObject::SaveToObj(const char* fn)
    {
        Array<int,3>* tlist = m_pTriangleList;
        edit();
        if(!IsLoaded())
        {
            printf("Error: Object is not initialized.\n",fn);
            return false;
        }
        
        FILE* fp = fopen(fn,"w");
        if(fp==NULL)
        {
            printf("Error: Failed opening %s to write\n",fn);
            return false;
        }
        
        fprintf(fp,"# %d vertices\n",m_newVertex);
        for(int i=0;i<m_newVertex;i++)
        {
            fprintf(fp,"v %f %f %f\n",  m_pVertexList[i].x,
                    m_pVertexList[i].y,
                    m_pVertexList[i].z);
            
        }
        
        fprintf(fp,"# %d triangles\n",m_newTriangle);
        for(int i=0;i<m_nTriangles;i++)
        {
            if ((tlist[i][0] == tlist[i][1])
                || (tlist[i][0] == tlist[i][2])
                || (tlist[i][1] == tlist[i][2]))//delete those contracted phase
                continue;
            fprintf(fp,"f %d %d %d\n",  m_pTriangleList[i][0] + 1,
                    m_pTriangleList[i][1] + 1,
                    m_pTriangleList[i][2] + 1);
        }
        
        fclose(fp);
        
        printf("Writing to %s successfully\n",fn);
        return true;
        
    }
}


