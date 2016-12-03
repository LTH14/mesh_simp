#pragma once
#include "Vec3f.h"
#include <vector>
#include <assert.h>
#include <set>
#include <list>

namespace SimpleOBJ
{
    
    
    template <typename T, int N> class Array
    {
    public:
        enum {_len = N};
        typedef T t_Val;
    public:
        T& operator[] (int i)
        {
            assert(i>=0&&i<N);
            return _p[i];
        }
        const T& operator[] (int i) const
        {
            assert(i>=0&&i<N);
            return _p[i];
        }
        
    protected:
        T _p[N];
    };
    
    struct Edge
    {
        int ID;
        float cost;
    };
    
    class EdgeSortCriterion {
    public:
        bool operator() (const Edge &a, const Edge &b) const {
            if(a.cost < b.cost)
                return true;
            else if(a.cost == b.cost)
            {
                if (a.ID < b.ID)
                    return true;
                else
                    return false;
            }
            else
                return false;
        }
    };

    
    class Qmetric
    {
    public:
        Qmetric()
        {
            for(int i = 0; i < 4; i++)
                for (int j = 0; j < 4; j++)
                    value[i][j] = 0;
        }
        ~Qmetric()
        {
            ;
        }
    public:
        float value[4][4];
    };
    
    
    class CSimpleObject
    {
    public:
        CSimpleObject(void);
        ~CSimpleObject(void);
        
    public:
        bool IsLoaded() { return m_pVertexList!=NULL;}
        
        void Destroy();
        bool LoadFromObj(const char* fn);
        bool SaveToObj(const char* fn);
        void contract();
        void compute_Q();
        void update_Q(int v1, int v2, Vec3f v1_value, Vec3f v2_value);
        void update_this_Q(int v1, int v2, Vec3f v1_value, Vec3f v2_value);
        float comp_err(int i, int j);
        Vec3f opt_target(int i, int j);
        int             m_nowTriangles;
        void init_adjlist();
        
    protected:
        bool Parse(FILE* fp);
        bool CheckParse(int nVertices,std::vector<Array<int,3> > & vecTriangles);
        
    protected:
        
        int             m_nVertices;
        int             m_nTriangles;
        
        Vec3f*          m_pVertexList;
        Array<int,3>*   m_pTriangleList;
        Qmetric*        m_pQList;
        std::set<Edge, EdgeSortCriterion> all_edges; // each element's i must be smaller than j
        int   id_num = 0;
    };
    
}

