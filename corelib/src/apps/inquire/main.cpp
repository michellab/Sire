
#include "qcustomplot.h"

#include <QApplication>
#include <QMainWindow>

#include "SireStream/streamdata.hpp"
#include "SireStream/metatype.h"

#include "SireBase/property.h"

#include "SireAnalysis/ti.h"

#include "SireError/errors.h"

#include <QDebug>

#include "tostring.h"

using namespace SireAnalysis;
using namespace SireStream;
using namespace SireBase;

int main(int argc, char **argv)
{
    QString file_to_load;

    if (argc > 0)
    {
        //load up the first s3 file in the argument list
        for (int i=1; i<argc; ++i)
        {
            QString arg(argv[i]);
            
            try
            {
                qDebug() << "Testing argument" << arg;
                FileHeader header = getDataHeader(arg);
                file_to_load = arg;
                qDebug() << "This is a valid s3 file :-)";
                qDebug() << header.requiredLibraries();
                qDebug() << header.dataTypes();
                break;
            }
            catch(...)
            {
                qDebug() << "Not a valid s3 file :-(";
            }
        }
    }

    QList<PropertyPtr> props;

    if (not file_to_load.isEmpty())
    {
        qDebug() << "Loading" << file_to_load;
        
        QList< boost::tuple<boost::shared_ptr<void>,QString> > objects;
        
        try
        {
            objects = load(file_to_load);
        }
        catch(const SireError::exception &e)
        {
            qDebug() << e.toString();
            return -1;
        }
        catch(...)
        {
            qDebug() << "Caught an unknown exception!";
            return -1;
        }
        
        qDebug() << "Loaded" << objects.count() << "objects";
        
        for (int i=0; i<objects.count(); ++i)
        {
            QString type_name = objects.at(i).get<1>();
            QString root_type = SireStream::registeredRoot(type_name);
            
            if (root_type == Property::typeName())
            {
                props.append( PropertyPtr(
                    static_cast<const Property*>(objects.at(i).get<0>().get())->clone() ) );
            }
        }
        
        objects.clear();
        
        foreach (PropertyPtr prop, props)
        {
            qDebug() << "Property" << prop.read().toString();
        }
    }
    
    QApplication *a = new QApplication(argc, argv);

    QCustomPlot *customPlot = new QCustomPlot();
    
    QVector<DataPoint> points;
    
    foreach (PropertyPtr prop, props)
    {
        if (prop.read().isA<TI>())
        {
            TI ti = prop.read().asA<TI>();
            qDebug() << ti.toString();
            qDebug() << ti[-1].toString();
            qDebug() << ti[-1].integrate().toString();
            points = ti.merge(ti.nIterations()/2, ti.nIterations()-1).integrate().values();
            qDebug() << Sire::toString(points);
            break;
        }
    }

    // create graph and assign data to it:
    customPlot->addGraph();
    
    QVector<double> x;
    QVector<double> y;
    
    double ymin = 1000000;
    double ymax = -1000000;
    
    foreach (DataPoint point, points)
    {
        x.append( point.x() );
        y.append( point.y() );
        
        ymin = qMin(ymin, point.y());
        ymax = qMax(ymax, point.y());
        
        qDebug() << point.x() << point.y();
    }
    
    double delta = (ymax - ymin);
    
    ymin -= 0.1*delta;
    ymax += 0.1*delta;
    
    qDebug() << ymin << ymax;
    
    customPlot->graph(0)->setData(x, y);
    // give the axes some labels:
    customPlot->xAxis->setLabel("lambda");
    customPlot->yAxis->setLabel("delta G");
    // set axes ranges, so we see all data:
    customPlot->xAxis->setRange(0, 1);
    customPlot->yAxis->setRange(ymin, ymax);
    customPlot->replot();

    customPlot->setWindowTitle("inquire");
    customPlot->resize(600,600);

    customPlot->show();

    return a->exec();
}
