package snetworktests;

import org.junit.runner.RunWith;
import org.junit.runners.Suite;
import org.junit.runners.Suite.SuiteClasses;

@RunWith(Suite.class)
@SuiteClasses({
    ConstantPopIntegratedTest.class,
    ConstantPopulationTest.class,
    NetworkParserTest.class,
    BirthHybridizationTest.class
})

public class AllTests {

}
