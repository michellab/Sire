
#include <QHash>

namespace Foo
{

    class Bar
    {
    public:
        Bar()
        {}

        ~Bar()
        {}

        bool operator==(const Bar &other) const
        {
            return true;
        }
    };

    /** The qHash function *must* be in the same namespace
        for the class for g++ 4.2 or else it won't be found.

        Try commenting out this function and then uncommenting
        out the qHash function below to see what I mean
    */
    uint qHash(const Bar &bar)
    {
        return 0;
    }
}

/** g++ 4.2 won't find this qHash function! (uncomment this function,
    then comment out the one above to see what I mean) */
//uint qHash(const Foo::Bar &bar)
//{
//    return 0;
//}

int main(int argc, const char **argv)
{
    QHash<Foo::Bar, Foo::Bar> foo_hash;

    foo_hash.insert( Foo::Bar(), Foo::Bar() );

    return 0;
}
