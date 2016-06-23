/*
 * ThreadUtils.hh
 *
 *  Created on: 16/01/2012
 *      Author: gdaviet
 */

#ifndef THREADUTILS_HH_
#define THREADUTILS_HH_

#include <boost/thread/thread.hpp>
#include <boost/thread/locks.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/condition_variable.hpp>

typedef boost::mutex MutexType ;
typedef boost::lock_guard< MutexType > LockGuard ;
typedef boost::unique_lock< MutexType > UniqueLock ;
typedef boost::condition_variable Condition ;

typedef boost::thread ThreadType ;

/*! Since boost::mutex are non copyable, this class allows other to have a member mutex
  without having to redefine their copy constructor and assignement operator
*/
class MutexWrapper
{
public:
    MutexWrapper()
    {}

    MutexWrapper( const MutexWrapper& )
    {}

    MutexWrapper& operator=( const MutexWrapper& )
    {
        return *this ;
    }

    MutexType& operator*()
    {
        return m_mutex ;
    }

private:
    MutexType m_mutex ;
};

class ThreadHandle
{
public:
    ThreadHandle()
    : m_running( false )
    {
    }

    ~ThreadHandle()
    {
        join() ;
    }

   template< typename DataT, void (DataT::*func)( void ) >
   void run( DataT* callee )
   {

        while ( m_running )
        {
            join() ;
        }

        LockGuard lock( m_mutex ) ;
        callable< DataT, func > f ;

        m_thread = boost::thread( f, callee ) ;
        m_running = true ;

    }

    template< typename DataT >
    void run( DataT* callee )
    {
        run< DataT, &DataT::operator() >( callee ) ;
    }



    void join()
    {
        LockGuard lock( m_mutex ) ;

        if( m_running )
        {
            m_thread.join() ;
            m_running = false ;
        }
    }

private:
    bool m_running ;
    boost::thread m_thread ;
    MutexType m_mutex ;

    template< typename DataT, void (DataT::*func)( void ) >
    struct callable
    {
        void operator()( DataT* callee )
        {
            (callee->*func)() ;
        }
    } ;


};

template<typename ContainerT, typename ReturnT, typename CallableT, typename ArgT>
void for_each( ContainerT& list, CallableT& obj, ReturnT(CallableT::*func)( ArgT ) )
{
    auto callee = std::bind1st( std::mem_fun( func ), &obj );
    for_each( list.begin(), list.end(), callee );
}

template<typename ContainerT, typename ReturnT, typename CallableT, typename ArgT>
void for_each( ContainerT& list, const CallableT& obj, ReturnT(CallableT::*func)( ArgT ) const )
{
    auto callee = std::bind1st( std::mem_fun( func ), &obj );
    for_each( list.begin(), list.end(), callee );
}

template<typename ContainerT, typename ReturnT, typename CallableT, typename ArgT>
void parfor( ContainerT& vec, CallableT& obj, ReturnT(CallableT::*func)( ArgT ) )
{
    const unsigned N = vec.size();
#pragma omp parallel for
    for ( unsigned i = 0; i < N; ++i )
    {
        ( obj.*func )( vec[i] );
    }
}

template<typename ContainerT, typename ReturnT, typename CallableT, typename ArgT>
void parfor( ContainerT& vec, const CallableT& obj, ReturnT(CallableT::*func)( ArgT ) const )
{
    const unsigned N = vec.size();
#pragma omp parallel for
    for ( unsigned i = 0; i < N; ++i )
    {
        ( obj.*func )( vec[i] );
    }
}

template<typename ContainerT, typename CallableT>
void parfor( ContainerT& vec, CallableT& obj )
{
    const unsigned N = vec.size();
#pragma omp parallel for
    for ( unsigned i = 0; i < N; ++i )
    {
        obj( vec[i] );
    }
}

#endif /* THREADUTILS_HH_ */
