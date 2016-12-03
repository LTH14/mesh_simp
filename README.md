网格简化步骤：
一个堆（set）存储所有边的cost。每个元素为二元组{int ID, float cost}
每一步：
1. pop堆直到pop出的ID有效。依据该ID找到相应的两个顶点v1,v2。
2. 根据面adjacency list，找到与v1, v2相邻的面，全部加入v1的面adjacency list。将这些面中的v2全部改为v1。
3. 将v1的值改为目标v的值。V2的值设为无效。
4. 根据边adjacency list，找到与v1，v2相邻的所有有效边ID。为其分配新的ID，将该新ID加入邻点的边adjacency list、v1的adjacency list，**计算该新边的cost，置入堆中。将原ID置为无效。
 
**
1. 对所有v1、v2的邻点，计算其新的Q（将邻点相邻的所有面再加一遍）。对v1同样计算其新的Q（此时面adjacency list已经完成，将v1相邻的所有面加一遍）。
2. 根据更新过的Q计算新边的cost。
 
所需的数据结构：
1. 一个堆（set）
2. 一个面adjacency list。横坐标为点，纵坐标为面。
3. 一个边adjacency list。横坐标为点，纵坐标为相邻边的ID
4. 一个边ID数组。横坐标为边ID，纵坐标为flag（是否有效），v1，v2（两个端点）。
5. 一个边ID个数的记录器
