package snetworktests;

import org.junit.runner.RunWith;
import org.junit.runners.Suite;
import org.junit.runners.Suite.SuiteClasses;

@RunWith(Suite.class)
@SuiteClasses({
    ConstantPopIOTest.class,
    ConstantPopulationTest.class,
    NetworkParserTest.class,
    YuleHybridModelTest.class
})

public class AllTests {

}
