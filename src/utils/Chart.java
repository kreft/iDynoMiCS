/**
 * \package utils
 * \brief Package of classes that perform utility functions in the process of running an iDynoMiCS Simulation
 * 
 * Package of classes that perform utility functions in the process of running an iDynoMiCS Simulation. This package is part of iDynoMiCS v1.2, governed by the 
 * CeCILL license under French law and abides by the rules of distribution of free software.  You can use, modify and/ or redistribute 
 * iDynoMiCS under the terms of the CeCILL license as circulated by CEA, CNRS and INRIA at the following URL  "http://www.cecill.info".
 */
package utils;


import java.io.File;
//import java.util.ArrayList;
//import java.util.Arrays;

import org.jfree.chart.ChartPanel;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;

import org.jfree.chart.axis.AxisLocation;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.CombinedDomainXYPlot;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.StandardXYItemRenderer;
import org.jfree.chart.renderer.xy.XYItemRenderer;
//import org.jfree.data.xy.XYDataset;
//import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.ui.ApplicationFrame;


@SuppressWarnings("serial")
/**
 * \brief Class used to represent simulation output on a graph. Assumed deprecated as never called in iDynoMiCS 1.2
 * 
 * Class used to represent simulation output on a graph. Assumed deprecated as never called in iDynoMiCS 1.2
 *
 */
public class Chart extends ApplicationFrame 
{

	/**
	 * Source of the data to be plotted
	 */
	private XYSeriesCollection[] _dataSource;
	
	/**
	 * Path to the file where this chart will be stored
	 */
	private String               _resultPath;
	
	/**
	 * The top-level chart object
	 */
	private CombinedDomainXYPlot parentPlot;
	
	/**
	 * The chart itself - making use of JFreeChart
	 */
	private JFreeChart           chart;
	
	/**
	 * Panel to be displayed alongside the chart
	 */
	private ChartPanel           panel;

	/**
	 * \brief Constructs the chart object on the screen, with title provided
	 * 
	 * Constructs the chart object on the screen, with title provided
	 * 
	 * @param title	Title to be used for this chart
	 */
	public Chart(final String title) 
	{

		super(title);

		parentPlot = new CombinedDomainXYPlot();

		parentPlot.setGap(10.0);
		parentPlot.setOrientation(PlotOrientation.VERTICAL);

		chart = new JFreeChart(title, JFreeChart.DEFAULT_TITLE_FONT, parentPlot, true);
		panel = new ChartPanel(chart, true, true, true, false, true);
		panel.setPreferredSize(new java.awt.Dimension(1024, 768));

		setContentPane(panel);
	}

	/**
	 * \brief Set the path to the location where this chart will be saved
	 * 
	 * Set the path to the location where this chart will be saved
	 * 
	 * @param resultPath	Path to a location on disk where this chart will be saved
	 */
	public void setPath(String resultPath) {
		_resultPath = resultPath;
	}

	/**
	 * \brief Initialise the data sets to be drawn in the created chart, and the chart axis
	 * 
	 * Initialise the data sets to be drawn in the created chart, and the chart axis
	 * 
	 * @param allGraphSets	Data to be plotted on the graph
	 * @param xLegend	X axis legend
	 * @param yLegend	Y axis legend
	 */
	public void init(XYSeriesCollection[] allGraphSets, String[] xLegend, String[] yLegend) {
		_dataSource = allGraphSets;
		XYPlot subplot;
		NumberAxis rangeAxisX, rangeAxisY;
		XYItemRenderer renderer;

		for (int iPlot = 0; iPlot<_dataSource.length; iPlot++) {
			renderer = new StandardXYItemRenderer();
			rangeAxisX = new NumberAxis(xLegend[iPlot]);
			rangeAxisY = new NumberAxis(yLegend[iPlot]);
			subplot = new XYPlot(allGraphSets[iPlot], rangeAxisX, rangeAxisY, renderer);
			subplot.setRangeAxisLocation(AxisLocation.BOTTOM_OR_LEFT);
			parentPlot.add(subplot, 1);
		}
	}

	/**
	 * \brief Update chart with the data provided in the input arguments
	 * 
	 * Update chart with the data provided in the input arguments
	 * 
	 * @param iChart	The chart being drawn
	 * @param iSeries	The data series on the chart
	 * @param x	The x point to add
	 * @param y	The y point to add
	 */
	public void updateChart(int iChart, int iSeries, double x, double y) {
		_dataSource[iChart].getSeries(iSeries).add(x, y);
	}

	/**
	 * \brief Update chart with the data provided in the input argument arrays
	 * 
	 * Update chart with the data provided in the input argument arrays
	 * 
	 * @param iChart	The chart being drawn
	 * @param iSeries	The data series on the chart
	 * @param x	Array of x points to add
	 * @param y	Array of y points to add
	 */
	public void updateChart(int iChart, int iSeries, double[] x, double[] y) {
		_dataSource[iChart].getSeries(iSeries).clear();
		for (int iPoint = 0; iPoint<x.length; iPoint++)
			_dataSource[iChart].getSeries(iSeries).add(x[iPoint], y[iPoint]);

	}

	/**
	 * \brief Refresh the graph on the screen and save to the result path location
	 * 
	 * Refresh the graph on the screen and save to the result path location
	 */
	public void repaintAndSave() 
	{
		this.repaint();
		try {
			ChartUtilities.saveChartAsPNG(new File(_resultPath+"solute.png"), chart, 1024, 768);
		} catch (Exception e) {
		}

	}
}
