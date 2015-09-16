#include <cstdio>
#include <vector>
#include <queue>
#include <climits>
#include <algorithm>

class Graph
{
    struct _Edge
    {
        unsigned int to;
        unsigned int capacity;
        int flow;
        unsigned int reverse;
        
        _Edge(int to, unsigned int capacity, unsigned int reverse):
            to(to),
            capacity(capacity),
            flow(0),
            reverse(reverse)
        {
        }
    };
    
    std::vector <std::vector <_Edge> > adjacencyLists;
    
    unsigned int bfs(unsigned int source, unsigned int sink)
    {
        std::vector<unsigned int> distance(adjacencyLists.size(), UINT_MAX);
        std::vector<unsigned int> parent(adjacencyLists.size(), source);
        parent[source] = UINT_MAX;
        distance[source] = 0;
        
        std::queue<unsigned int> q;
        q.push(source);
        
        while (q.size()) 
        {
            unsigned int v = q.front();
            q.pop();
            if (v == sink)
                break;
            for (auto const &e: adjacencyLists[v])
            {
                if ((e.capacity > e.flow) && (distance[e.to] > distance[v] + 1))
                {
                    distance[e.to] = distance[v] + 1;
                    parent[e.to] = v;
                    q.push(e.to);
                }
            }
        }
        
        if (distance[sink] == UINT_MAX)
        {
            return 0;
        }
        
        std::vector <unsigned int> way;
        
        way.push_back(sink);
        
        unsigned int minRemainingCapacity = UINT_MAX;
        
        
        while (parent[sink] != UINT_MAX)
        {
            unsigned int newSink = parent[sink];            
            
            unsigned int currentRemainingCapacity = 0;
            for (auto const &e: adjacencyLists[newSink])
            {
                if (e.to == sink)
                {
                    currentRemainingCapacity += e.capacity - e.flow;
                }
            }
            minRemainingCapacity = std::min(minRemainingCapacity, currentRemainingCapacity);
            
            sink = newSink;
            way.push_back(sink);
        }
        
        std::reverse(way.begin(), way.end());
        
        for (unsigned int i = 0; i + 1 < way.size(); ++i)
        {
            unsigned int currentFlow = minRemainingCapacity;
            
            for (auto &e: adjacencyLists[way[i]])
            {
                if (e.to == way[i + 1])
                {
                    unsigned int currentPush = std::min(e.capacity - e.flow, currentFlow);
                    e.flow += currentPush;
                    adjacencyLists[e.to][e.reverse].flow -= currentPush;
                    currentFlow -= currentPush;
                }
            }
        }
        
        return minRemainingCapacity;
    }
    
public:
    struct Edge
    {
        unsigned int from;
        unsigned int to;
        unsigned int capacity;
        
        Edge(unsigned int from, unsigned int to, unsigned int capacity):
            from(from),
            to(to),
            capacity(capacity)
        {
        }
        
        Edge()
        {
        }
    };
    
    Graph(unsigned int verticesCount, const std::vector <Edge> &edgeList)
    {
        adjacencyLists.resize(verticesCount);
        for (auto const &e: edgeList)
        {
            adjacencyLists[e.from].push_back(_Edge(e.to, e.capacity, adjacencyLists[e.to].size()));
            adjacencyLists[e.to].push_back(_Edge(e.from, 0, adjacencyLists[e.from].size() - 1));
        }
    }
    
    unsigned long long maxFlow(unsigned int source, unsigned int sink)
    {
        if (source == sink)
        {
            return ULLONG_MAX;
        }
        
        unsigned long long flowValue = 0;
        unsigned int deltaFlow;
        while ((deltaFlow = bfs(source, sink)))
        {
            flowValue += deltaFlow;
        }
        return flowValue;
    }
};


int main()
{
    unsigned int n, m;
    scanf("%u %u", &n, &m);
    
    std::vector<Graph::Edge> edges(m);
    
    for (unsigned int i = 0; i < m; ++i)
    {
        scanf("%u %u %u", &edges[i].from, &edges[i].to, &edges[i].capacity);
        --edges[i].from, --edges[i].to;
    }
    
    Graph g(n, edges);
    
    printf("%llu\n", g.maxFlow(0, n - 1));
    return 0;
}