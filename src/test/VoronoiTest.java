package test;

import java.util.Arrays;
import java.util.List;
import java.util.ListIterator;

public class VoronoiTest
{
	public static void main(String[] args) 
	{
		
		List<String> numbers = Arrays.asList("zero", "one", "two");
		ListIterator<String> it = numbers.listIterator();
		String[] out = new String[numbers.size()];
		while (it.hasNext())
		{
			out[it.nextIndex()] = it.next();
			//System.out.println(it.nextIndex() + " " + it.next());
			if ( out[it.previousIndex()].contains("n") )
				out[it.previousIndex()] = "blah";
			//System.out.println(it.previousIndex() + " " + it.previous());
			//System.out.println(it.nextIndex() + " " + it.next());
			//System.out.println(" --- ");
		}
		for ( String s : out )
			System.out.println(s);
		
	}
}
